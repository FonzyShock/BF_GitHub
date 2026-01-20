% Optimal_UAV_Path_Casadi.m
% Main optimization script for UAV Trajectory Generation using 
% Mixed Bernstein-Fourier Parameterization.
addpath('C:\Program Files\MATLAB\casadi-3.7.1-windows64-matlab2018b')
rng(123);                                               % seed random number
data = struct();                                        % Initialize the struct to hold all data

% -------------------------------------------------------------------------
% 1. Load Input Arguments or Set Defaults
% -------------------------------------------------------------------------
if exist('inputArgs', 'var')
    fn = fieldnames(inputArgs);
    for i = 1:numel(fn)
        assignin('base', fn{i}, inputArgs.(fn{i}));
    end
    
    % Camera struct (res_x,res_y,total_MP)
    CAM     = camera_resolution_catalog();
    cam_idx = max(1, min(numel(CAM), round(inputArgs.camera_idx)));
    camera  = CAM(cam_idx); 
    res_x   = camera.res_x;
    res_y   = camera.res_y;
    isDOE   = true;  
else
    clc; clear; close all;
    
    % Default values
    num_UAVs            = 1;                                
    boresight           = 45;                               
    res_x               = 2592;
    res_y               = 1944;
    isDOE               = false;                            
end

% Field of View (FOV) Settings
el_lim_up   = pi/2;
el_lim_dn   = pi/2;
az_lim      = pi;
sample_rate = 10;                               % Hz
V_A         = 25;                               % m/s
control_limit = 60;                             % Degrees
tf          = 30;                               % Seconds

% -------------------------------------------------------------------------
% 2. Problem Setup & Bernstein-Fourier Basis Definition
% -------------------------------------------------------------------------
data.Nx          = 3;                           % States: [x, y, psi]
data.Nu          = 1;                           % Control: [tan(phi)]
data.tf          = tf;                          
data.alpha_min   = -0.9;                        % Time scaling bounds
data.alpha_max   = 1.3;                         
data.fs          = sample_rate;                 

% --- Basis Parameters ---
data.N_B         = 30;                          % Bernstein Order (Polynomial/Trend)
data.N_F         = 12;                          % Fourier Order (Trigonometric/Orbit)

% TOTAL PARAMETERS: (Bernstein+1) + (2*Fourier)
% Note: We exclude the DC component from Fourier as Bernstein handles the constant offset.
data.N_tot       = (data.N_B + 1) + (2 * data.N_F); 

data.Nt          = data.tf * sample_rate + 1;           % Points for Plotting
data.time_opt    = linspace(0, 1, data.N_tot);          % Nodes for Optimization (Collocation)
data.time        = linspace(0, 1, data.Nt);             % Nodes for Final Evaluation
data.sample_grid = (0:data.N_tot)/data.N_tot;           
data.tol         = 1e-10;                               

% -------------------------------------------------------------------------
% 3. Generate Basis Matrices (Pre-computed)
% -------------------------------------------------------------------------
% Optimization Grid (used in constraints)
data.BN          = bernsteinMatrix(data.N_B, data.time_opt);
data.FN          = FourierMatrix(data.N_F, data.time_opt); 
data.DBm         = bernsteinMatrixDerivative(data.N_B, data.time_opt);
data.DFm         = FourierDiff(data.N_F, data.time_opt);

% Evaluation Grid (used in plotting/analysis)
data.BN_eval     = bernsteinMatrix(data.N_B, data.time);
data.FN_eval     = FourierMatrix(data.N_F, data.time);
data.DBm_eval    = bernsteinMatrixDerivative(data.N_B, data.time);
data.DFm_eval    = FourierDiff(data.N_F, data.time);

% -------------------------------------------------------------------------
% 4. Physics & Constraints
% -------------------------------------------------------------------------
data.lambda           = 1/num_UAVs;             
data.g                = 9.81;                   
data.rCom             = 15;                     
data.r_NFZ            = 100;                    
data.V_A              = V_A;                    
data.sample_rate      = sample_rate;            
data.dps              = data.V_A / data.sample_rate;         
bearing_error         = deg2rad(2);             
data.bearing_error    = bearing_error;          
sigma                 = bearing_error^2;        
data.Sigma            = diag(repmat(sigma,1,4));
data.Q                = 0.7 * eye(3);           
data.az_lim           = az_lim * pi / 180;      
data.el_lim_up        = el_lim_up * pi / 180;   
data.el_lim_dn        = el_lim_dn * pi /180;    
data.control_limit    = deg2rad(control_limit); 
data.MissionAreaBound = [-50000; 50000; -50000; 50000]; 

% Gimbal Settings
data.isGimballed      = false;                  
data.u_psi_g_max      = deg2rad(300);           
data.u_psi_g_min      = -data.u_psi_g_max;      
data.u_theta_g_max    = deg2rad(300);           
data.u_theta_g_min    = -data.u_theta_g_max;    

if data.isGimballed
    data.Nx           = 7; % [x, y, psi, cPsi_g, sPsi_g, cTheta_g, sTheta_g]
    data.Nu           = 3; % [tan_phi, u_psi_g, u_theta_g]
    num_UAVs          = 1;
end

% -------------------------------------------------------------------------
% 5. Mission Setup (Targets & Initial Conditions)
% -------------------------------------------------------------------------
data.num_targets      = 1;                      
data.priority         = 0.9;                    
data.time_weight      = 4e-1;                   

% Target Trajectory
x0_T    = [5 0] * 1e2;                          
data.V_T = 0;                                   
data.theta_T = 45;                              
data.cueing.target = 50;                        
data.swell_height = 1;                          

x_T = x0_T(1) + data.V_T * cosd(data.theta_T) * data.time_opt * data.tf;
y_T = x0_T(2) + data.V_T * sind(data.theta_T) * data.time_opt * data.tf;
z_T = zeros(1, data.N_tot);
data.x_T = [x_T; y_T; z_T];

var_xy0 = data.cueing.target;
var_z0  = data.swell_height;  
data.P0_AT = diag([var_xy0^2, var_xy0^2, var_z0^2]);

% Boat Trajectory
x0_B = [0 0] * 1e2;                             
data.V_B = 0;                                   
data.theta_B = 30;                              
data.cueing.boat = 10;                          

x_B = x0_B(1) + data.V_B * cosd(data.theta_B) * data.time_opt * data.tf;
y_B = x0_B(2) + data.V_B * sind(data.theta_B) * data.time_opt * data.tf;
z_B = zeros(1, data.N_tot);
data.X_B = [x_B; y_B; z_B];

var_xy0 = data.cueing.boat;
data.P0_AB = diag([var_xy0^2, var_xy0^2, var_z0^2]);

% UAV Initial Conditions
x0      = x0_B(1);
y0      = x0_B(2);
x_range = [x0 - 1, x0 + 1];
y_range = [y0 - 1, y0 + 1];
x_init  = x_range(1) + (x_range(2) - x_range(1)) * rand(num_UAVs, 1);
y_init  = y_range(1) + (y_range(2) - y_range(1)) * rand(num_UAVs, 1);
y_init  = sort(y_init, 1, "descend");            

UAV = zeros(num_UAVs, 9);
for idx = 1:num_UAVs
    Altitude      = 90 + idx*10; 
    x_0           = x_init(idx);
    y_0           = y_init(idx);
    dx_target     = x0_T(1) - x_0;
    dy_target     = x0_T(2) - y_0;
    target_angle  = atan2(dy_target, dx_target);
    
    launch_window = 45;   
    psi_init      = deg2rad(launch_window) * (rand - 0.5);
    Psi_0         = target_angle + psi_init;
    
    % Ensure forward flight
    u_hat         = [cos(Psi_0); sin(Psi_0)];
    d_vec         = [dx_target; dy_target];
    if dot(u_hat, d_vec) < 0, Psi_0 = Psi_0 + pi; end
    Psi_0         = wrapToPi(Psi_0);
    Phi_0         = 0;

    % Final Condition
    x_f           = x_B(end) + 0.5*data.rCom; 
    y_f           = y_B(end) + (0.5 + 0.5 * rand()) * data.rCom;
    dx_final      = x_B(end) - x_f;
    dy_final      = y_B(end) - y_f;
    Psi_f         = atan2(dy_final, dx_final);              
    Phi_f         = 0;

    UAV(idx,:)    = [Altitude, x_0, y_0, Psi_0, Phi_0, x_f, y_f, Psi_f, Phi_f];
end

data.num_uav = size(UAV, 1);
for uavIdx = 1:data.num_uav
    idx = num2str(uavIdx);
    data.(['z_A' idx])  = UAV(uavIdx, 1);                                   
    data.(['x0_A' idx]) = [UAV(uavIdx, 2), UAV(uavIdx, 3), UAV(uavIdx, 4)]; 
    data.(['u0' idx])   = tan(UAV(uavIdx, 5));                              
    data.(['xf_A' idx]) = [UAV(uavIdx, 6), UAV(uavIdx, 7), UAV(uavIdx, 8)]; 
    data.(['uf' idx])   = tan(UAV(uavIdx, 9));                              
end

% -------------------------------------------------------------------------
% 6. Optimization Execution
% -------------------------------------------------------------------------
import casadi.* 
opti = casadi.Opti(); 

opts = struct();      
opts.ipopt.acceptable_iter  = 1e2;           
opts.ipopt.max_iter         = 5e3;           
opts.ipopt.acceptable_tol   = 1e-3;          
opts.ipopt.tol              = 1e-4;          
opts.ipopt.constr_viol_tol  = 1e-6;          
opts.ipopt.mu_strategy      = 'adaptive';    
opts.ipopt.print_level      = 5;             
opts.print_time             = false;
data.opts                   = opts;

data.warmstart.clearCache   = false;
data.warmstart.useCache     = false;
data.isDOE                  = isDOE;
data.plotinit               = true;
data.plotfig                = true; 
data.reverse_flight_dir     = false; 
data_nonCoop                = data;

% --- A. Solve Non-Cooperative ---
for uavIdx = 1:data.num_uav
    idx     = num2str(uavIdx);
    xout    = ['xOut' idx];
    jout    = ['Jout' idx];
    
    data_nonCoop.idx = idx; 
    [results_nonCoop.(xout), results_nonCoop.(jout), data_nonCoop] = ...
        nonCooperativeUAVs(data_nonCoop, uavIdx);
end

[tr_crlb_AB, tr_crlb_AT, data_nonCoop] = calc_CRLB_noncoop(data_nonCoop, results_nonCoop);

% --- B. Solve Cooperative ---
data_Coop = data_nonCoop;
opts.ipopt.acceptable_tol   = 1e-4;          
opts.ipopt.tol              = 1e-5;          
opts.ipopt.constr_viol_tol  = 1e-5;          
data_Coop.opts              = opts;

[results_Coop.xout, results_Coop.jout, data_Coop] = cooperativeUAVs(data_Coop, results_nonCoop); 
[tr_crlb_AB, tr_crlb_AT, data_Coop] = calc_CRLB_coop(data_Coop, results_Coop);

tr_CRLB_Coop = struct('tr_crlb_AB', tr_crlb_AB, 'tr_crlb_AT', tr_crlb_AT);
converged    = struct('convergeFlag', data_Coop.convergeFlag);

% -------------------------------------------------------------------------
% 7. Plotting & EKF
% -------------------------------------------------------------------------
if data.plotfig && isequal(isDOE, false)
    plotResults(data_nonCoop, results_nonCoop, data_Coop, results_Coop);
end

% try
%     run('Federated_EKF.m');
%     DOE_Input.rmse = RMSE_MC;
% catch
%     fprintf('Federated_EKF.m not found on path.\n');
% end


% =========================================================================
%  BASIS FUNCTIONS (MUST Match initPath.m)
% =========================================================================

function BN = bernsteinMatrix(N, time)
    t = (time - time(1)) / (time(end) - time(1));
    t = t(:);
    BN = zeros(length(t), N + 1);
    
    for k = 0:N
        log_coeff = gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1);
        % Safe computation for t=0 and t=1
        term = exp(log_coeff + k .* log(t + 1e-16) + (N - k) .* log(1 - t + 1e-16));
        BN(:, k + 1) = term;
    end
    % Exact endpoints
    BN(1, 1) = 1; BN(1, 2:end) = 0;
    BN(end, end) = 1; BN(end, 1:end-1) = 0;
end

function BN_dot = bernsteinMatrixDerivative(N, time)
    t = (time - time(1)) / (time(end) - time(1));
    t = t(:);
    BN_minus1 = bernsteinMatrix(N-1, t);
    
    BN_left  = [zeros(length(t), 1), BN_minus1];
    BN_right = [BN_minus1, zeros(length(t), 1)];
    
    BN_dot = N * (BN_left - BN_right);
end

function FN = FourierMatrix(N, time)
    % Fourier Basis without DC term
    % Returns [cos(wt), sin(wt), cos(2wt), sin(2wt)...]
    t = time(:); 
    Nt = numel(t); 
    T = time(end);
    if T==0, T=1; end
    if max(time) <= 1.0001, T = 1; end
    
    w_base = 2*pi*(t/T);
    FN = zeros(Nt, 2*N);
    
    for k = 1:N
        w = w_base * k;
        col_cos = (k-1)*2 + 1;
        col_sin = (k-1)*2 + 2;
        FN(:, col_cos) = cos(w);
        FN(:, col_sin) = sin(w);
    end
end

function Dm = FourierDiff(N, time)
    % Derivative of Fourier Basis
    t = time(:); 
    Nt = numel(t); 
    T = time(end);
    if T==0, T=1; end
    if max(time) <= 1.0001, T = 1; end
    
    w_base = 2*pi*(t/T);
    Dm = zeros(Nt, 2*N);
    
    for k = 1:N
        w = w_base * k;
        rate = k * (2*pi/T);
        
        col_cos = (k-1)*2 + 1;
        col_sin = (k-1)*2 + 2;
        
        % d/dt(cos) = -rate*sin
        % d/dt(sin) =  rate*cos
        Dm(:, col_cos) = -rate .* sin(w);
        Dm(:, col_sin) =  rate .* cos(w);
    end
end