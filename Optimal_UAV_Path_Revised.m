% Optimal_UAV_Path_Revised.m
% =========================================================================
% Optimal Observer Trajectory Generation for Moving Target Localization
% Using Mixed Bernstein-Fourier Approximants (3D + FOV + Turn Rate)
%
% Adapted from moving_target_Liraz.m to support:
%   - 3D Physics & Turn-rate control (Dubins-like)
%   - Field of View (FOV) constraints
%   - Boat + Target FIM fusion
%
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  COMMON PROBLEM PARAMETERS
%  ========================================================================
fprintf('=========================================================================\n');
fprintf('  3D MOVING TARGET OBSERVER TRAJECTORY OPTIMIZATION\n');
fprintf('=========================================================================\n\n');

% --- 1. Load Defaults (from Optimal_UAV_Path_Casadi.m) ---
num_UAVs    = 1;
tf          = 50;                       % s, scenario duration
sample_rate = 10;                       % Hz
V_A         = 25;                       % m/s, UAV constant velocity
g           = 9.81;
control_lim = deg2rad(60);              % Max bank angle
z_A         = 100;                      % m, UAV altitude

% Map Boundaries
MissionAreaBound = [-50000; 50000; -50000; 50000]; 

% Targets (Boat + Target)
x0_T    = [5 0] * 1e2;
x0_B    = [0 0] * 1e2;
V_T     = 0;  theta_T = 45;
V_B     = 0;  theta_B = 30;

% Measurement Parameters
Sigma   = diag([deg2rad(2)^2, deg2rad(2)^2]); % [az^2, el^2]
Q       = 0.7 * eye(3);
priority= 0.9;                          % Weight for Target vs Boat

% Constraint Parameters
r_NFZ   = 100;                          % m, no-fly zone radius
rCom    = 15;                           % m, comms range (not used in barrier, but good to know)
kappa   = 1e-2;                         % Barrier parameter

% FOV Parameters
az_lim      = deg2rad(100);             % Approx from catalog
el_lim_up   = deg2rad(45);
el_lim_dn   = deg2rad(45);

% Approximation orders
N_B = 30;                               % Bernstein polynomial order
N_F = 12;                               % Fourier series order

% State and Control dimensions
% States: [x, y, psi], Controls: [tan(phi)]
Nx = 3; 
Nu = 1; 

%% ========================================================================
%  PART 1: MIXED BERNSTEIN-FOURIER OPTIMIZATION
%  ========================================================================
fprintf('------------------------------------------------------------------\n');
fprintf('  PART 1: Mixed Bernstein-Fourier Optimization\n');
fprintf('------------------------------------------------------------------\n');

% Total number of coefficients
N_tot = N_B + 1 + 2*N_F; 

% Time nodes for Optimization
Nt_BF   = N_B + 1; % Or use N_tot for collocation
time_BF = linspace(0, 1, N_tot); % Collocation points

% Construct basis matrices
BN_BF  = bernsteinMatrix(N_B, time_BF);
FN_BF  = FourierMatrix(N_F, time_BF);
DBm_BF = bernsteinMatrixDerivative(N_B, time_BF);
DFm_BF = FourierDiff(N_F, time_BF);

fprintf('Basis matrix condition number: %.2e\n', cond([BN_BF, FN_BF]));

% --- Initial Guess Generation ---
% We use a simple circular orbit guess around the midpoint or target
% to avoid needing external initPath dependencies, ensuring "standalone" robustness.
fprintf('Generating Initial Guess...\n');

Cx_BF = zeros(N_tot, Nx);
Cu_BF = zeros(N_tot, Nu);

% Simple strategy: Fly straight to target area, then circle
t_sec = time_BF * tf;
start_pos = [x0_B(1), x0_B(2)]; % Start near boat
end_pos   = [x0_T(1), x0_T(2)]; % End near target

for k = 1:N_tot
    alpha = time_BF(k);
    % Linear interp position
    pos_k = (1-alpha)*start_pos + alpha*end_pos;
    % Add some "orbit" wiggle to make it feasible-ish
    orbit_radius = 200;
    pos_k = pos_k + [orbit_radius*cos(2*pi*alpha), orbit_radius*sin(2*pi*alpha)];
    
    % Approx heading (tangent to circle)
    psi_k = atan2(end_pos(2)-start_pos(2), end_pos(1)-start_pos(1)) + pi/2 + 2*pi*alpha;
    
    Cx_BF(k, :) = [pos_k(1), pos_k(2), unwrap(psi_k)];
    Cu_BF(k, :) = 0; % Zero bank angle initially
end

% Reshape for optimizer: X = [Cx_vec; Cu_vec; alpha_time]
% Note: Liraz's code puts states then controls. 
% We add alpha (time scaling) as the last var to match your original Casadi logic if needed,
% but Liraz's code fixes tf. We will FIX tf here to match Liraz's structure "almost exactly".
X0_BF = [reshape(Cx_BF, [Nx * N_tot, 1]); reshape(Cu_BF, [Nu * N_tot, 1])];

% --- Constraints Setup ---
% Initial Conditions (Equality)
% State at t=0 must match x0
x0_A = [x0_B(1)-100, x0_B(2), 0]; % Offset slightly
u0   = 0;

% Matrix to extract t=0 state from coefficients
Basis_0 = [BN_BF(1,:), FN_BF(1,:)]; 
Aeq_BF  = zeros(Nx, (Nx+Nu)*N_tot);
beq_BF  = zeros(Nx, 1);

for i = 1:Nx
    col_start = (i-1)*N_tot + 1;
    col_end   = i*N_tot;
    Aeq_BF(i, col_start:col_end) = Basis_0;
    beq_BF(i) = x0_A(i);
end

% Bounds
lb_BF = -inf(size(X0_BF));
ub_BF =  inf(size(X0_BF));

% Enforce control limits on coefficients (convex hull property holds for B-splines, 
% less so for Fourier, but serves as a good box constraint)
u_idx_start = Nx * N_tot + 1;
lb_BF(u_idx_start:end) = -tan(control_lim);
ub_BF(u_idx_start:end) =  tan(control_lim);

% Optimizer options (matching Liraz + Fmincon)
options_BF = optimoptions(@fmincon, ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 2e5, ...
    'OptimalityTolerance', 1e-4,...
    'ConstraintTolerance', 1e-4, ...
    'StepTolerance', 1e-8, ...
    'MaxIterations', 3000, ...
    'Display', 'iter');

% --- Structs for Target/Boat Data (for Cost Func) ---
missionData.x0_T    = x0_T;
missionData.V_T     = V_T; missionData.theta_T = theta_T;
missionData.x0_B    = x0_B;
missionData.V_B     = V_B; missionData.theta_B = theta_B;
missionData.z_A     = z_A;
missionData.Sigma   = Sigma;
missionData.az_lim  = az_lim;
missionData.el_lim_up = el_lim_up;
missionData.el_lim_dn = el_lim_dn;
missionData.priority  = priority;

% Define Cost and Constraints
costFunc_BF = @(X) costFunc_3D_FOV(X, N_B, N_F, N_tot, Nx, Nu, ...
    missionData, tf, r_NFZ, kappa, BN_BF, FN_BF, time_BF);

nonlcon_BF = @(X) nonlcon_3D_Dynamics(X, N_B, N_F, N_tot, Nx, Nu, ...
    V_A, g, tf, DBm_BF, DFm_BF, BN_BF, FN_BF);

% Run Optimization
fprintf('\nStarting Mixed BF optimization...\n');
tic;
[xOut_BF, Jout_BF, exitflag_BF, output_BF] = ...
    fmincon(costFunc_BF, X0_BF, [], [], Aeq_BF, beq_BF, lb_BF, ub_BF, nonlcon_BF, options_BF);
time_BF = toc;

fprintf('Mixed BF completed in %.2f seconds\n', time_BF);
fprintf('  Cost (CRLB + barrier): %.5e\n', Jout_BF);
fprintf('  Exit flag: %d\n', exitflag_BF);
fprintf('  Iterations: %d\n', output_BF.iterations);

% Store Results
results.mixedBF.xOut = xOut_BF;
results.mixedBF.time = time_BF;
results.mixedBF.exitflag = exitflag_BF;

%% ========================================================================
%  PART 2: EVALUATE AND PLOT (Based on Part 3/6 of Liraz's code)
%  ========================================================================
fprintf('\n------------------------------------------------------------------\n');
fprintf('  PART 2: Evaluate and Generate Figures\n');
fprintf('------------------------------------------------------------------\n');

% Evaluation Grid (Finer)
Nt_eval = 200;
time_eval = linspace(0, 1, Nt_eval);
time_eval_sec = time_eval * tf;

BN_eval = bernsteinMatrix(N_B, time_eval);
FN_eval = FourierMatrix(N_F, time_eval);

[X_A, Y_A, Psi_A, tan_Phi] = getStates_3D(xOut_BF, N_B, N_F, N_tot, Nx, Nu, BN_eval, FN_eval);

% Reconstruct Target/Boat Trajectories
X_T = x0_T(1) + V_T*cosd(theta_T)*time_eval_sec;
Y_T = x0_T(2) + V_T*sind(theta_T)*time_eval_sec;
X_B = x0_B(1) + V_B*cosd(theta_B)*time_eval_sec;
Y_B = x0_B(2) + V_B*sind(theta_B)*time_eval_sec;

% Figure 1: Top-Down Trajectory
figure('Name', 'Optimal Trajectory 3D', 'Position', [100, 100, 800, 600]);
hold on; grid on; axis equal;

% Plot UAV
scatter(X_A, Y_A, 20, time_eval_sec, 'filled'); 
colorbar; ylabel(colorbar, 'Time (s)');
plot(X_A, Y_A, 'b-', 'LineWidth', 1.5);

% Plot Targets
plot(X_T, Y_T, 'r--', 'LineWidth', 2, 'DisplayName', 'Target');
plot(X_B, Y_B, 'g--', 'LineWidth', 2, 'DisplayName', 'Boat');

% Plot NFZ
theta = linspace(0, 2*pi, 100);
plot(x0_T(1) + r_NFZ*cos(theta), x0_T(2) + r_NFZ*sin(theta), 'r:', 'LineWidth', 1.5, 'DisplayName', 'NFZ');

xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('Optimal UAV Trajectory (Cost: %.4f)', Jout_BF));
legend('Location','best');

% Figure 2: States
figure('Name', 'States vs Time', 'Position', [150, 150, 800, 600]);
subplot(3,1,1); plot(time_eval_sec, X_A); ylabel('X [m]'); grid on; title('States');
subplot(3,1,2); plot(time_eval_sec, Y_A); ylabel('Y [m]'); grid on;
subplot(3,1,3); plot(time_eval_sec, rad2deg(Psi_A)); ylabel('Psi [deg]'); xlabel('Time [s]'); grid on;

% Figure 3: Controls
figure('Name', 'Controls', 'Position', [200, 200, 800, 300]);
plot(time_eval_sec, rad2deg(atan(tan_Phi)), 'LineWidth', 2);
yline(rad2deg(control_lim), 'r--'); yline(-rad2deg(control_lim), 'r--');
ylabel('Bank Angle [deg]'); xlabel('Time [s]'); grid on; title('Control Input');


%% ========================================================================
%  HELPER FUNCTIONS (Embedded for completeness)
%  ========================================================================

%% Cost Function: 3D FIM Trace + FOV + Barrier
function J = costFunc_3D_FOV(X, N_B, N_F, N_tot, Nx, Nu, data, tf, r_NFZ, kappa, BN, FN, time)
    
    % 1. Extract States on optimization grid
    [X_A, Y_A, Psi_A, ~] = getStates_3D(X, N_B, N_F, N_tot, Nx, Nu, BN, FN);
    Z_A = data.z_A; 
    
    % 2. Target & Boat Trajectories (Vectorized over time)
    t_sec = time(:) * tf;
    Nt    = length(t_sec);
    
    X_T = data.x0_T(1) + data.V_T * cosd(data.theta_T) * t_sec;
    Y_T = data.x0_T(2) + data.V_T * sind(data.theta_T) * t_sec;
    
    X_B = data.x0_B(1) + data.V_B * cosd(data.theta_B) * t_sec;
    Y_B = data.x0_B(2) + data.V_B * sind(data.theta_B) * t_sec;
    
    % 3. Calculate FIM for Target (Vectorized 3D)
    % Relative position
    dx = X_T - X_A; 
    dy = Y_T - Y_A; 
    dz = -Z_A; % Assuming target at z=0
    
    [FIM_trace_T, dH_T] = compute_vectorized_FIM_trace(dx, dy, dz, Psi_A, data);
    
    % 4. Calculate FIM for Boat
    dx_b = X_B - X_A; 
    dy_b = Y_B - Y_A; 
    
    [FIM_trace_B, ~] = compute_vectorized_FIM_trace(dx_b, dy_b, dz, Psi_A, data);
    
    % 5. Barrier for NFZ (Target only)
    % Penalize if dist < r_NFZ. Liraz uses max(0, r - dH)^2
    barrier_term = kappa * norm(max(0, r_NFZ - dH_T))^2;
    
    % 6. Total Cost
    % Weighted sum of traces (Inverse FIM is roughly proportional to 1/trace(FIM) for optimization purposes, 
    % but strictly we sum FIMs then invert. Here we minimize trace(inv(FIM)) approximated or sum of individual traces)
    % Liraz minimizes trace(inv(FIM)). 
    
    % Since we have a sum over time, we sum the cost at each step
    J = data.priority * sum(FIM_trace_T) + (1-data.priority) * sum(FIM_trace_B) + barrier_term;
end

function [J_trace, dH] = compute_vectorized_FIM_trace(dx, dy, dz, Psi, data)
    % Vectorized FIM calculation with Sigmoid FOV
    
    rxy2 = dx.^2 + dy.^2 + 1e-6;
    rxy  = sqrt(rxy2);
    r2   = rxy2 + dz^2;
    dH   = sqrt(r2);
    
    % --- FOV Check (Vectorized InFOV logic) ---
    % 1. Body Frame Transformation
    cPsi = cos(Psi); sPsi = sin(Psi);
    % x_body = cPsi.*dx + sPsi.*dy;
    % y_body = -sPsi.*dx + cPsi.*dy;
    % For Azimuth check: atan2(y_body, x_body)
    % Actually InFOV uses rel_az in body frame
    
    x_b =  cPsi.*dx + sPsi.*dy;
    y_b = -sPsi.*dx + cPsi.*dy;
    z_b = dz; % Flat flight
    
    rel_az = atan2(y_b, x_b);
    rel_el = atan2(-z_b, sqrt(x_b.^2 + y_b.^2));
    
    % Sigmoid Soft FOV
    k_sig = 10;
    sig = @(x, thr) 1./(1 + exp(k_sig*(x - thr))); % Smooth step down
    
    % Azimuth check (abs(az) < lim) -> (az < lim) & (az > -lim)
    in_az = 1./(1 + exp(k_sig*(abs(rel_az) - data.az_lim/2))); 
    
    % Elevation check
    in_el = sig(rel_el, data.el_lim_up) .* (1 - sig(rel_el, data.el_lim_dn));
    
    % Visibility Weight
    w_vis = in_az .* in_el;
    
    % --- FIM Elements (3D Bearing Only) ---
    % Jacobian of [az; el] wrt [x,y,z]
    % H_az = [-dy/rxy2, dx/rxy2, 0]
    % H_el = [dx*dz/(r2*rxy), dy*dz/(r2*rxy), -rxy/r2]
    
    % We compute trace(inv(FIM)). For 3D position [x,y,z], FIM is 3x3.
    % To match Liraz's speed, we do this analytically or simplified.
    % Simplification: We minimize 1/det(FIM) or sum(1./eigenvalues).
    % Alternatively, define Cost = 1 / (Information_Gain).
    
    % Let's use the explicit structure:
    % FIM = Sum_k ( H_k' * R^-1 * H_k * w_vis_k )
    % This is expensive to invert inside loop.
    % Fallback: Maximize determinants/traces of *instantaneous* FIMs? 
    % Liraz sums the FIM elements over time, THEN inverts.
    
    % Construct Cumulative FIM Elements
    inv_sig_az = 1/data.Sigma(1,1);
    inv_sig_el = 1/data.Sigma(2,2);
    
    % Element-wise terms
    h_az_1 = -dy./rxy2; h_az_2 = dx./rxy2;
    h_el_1 = (dx.*dz)./(r2.*rxy); h_el_2 = (dy.*dz)./(r2.*rxy); h_el_3 = -rxy./r2;
    
    % FIM = [a b d; b c e; d e f]
    a = w_vis .* (inv_sig_az * h_az_1.^2 + inv_sig_el * h_el_1.^2);
    b = w_vis .* (inv_sig_az * h_az_1.*h_az_2 + inv_sig_el * h_el_1.*h_el_2);
    c = w_vis .* (inv_sig_az * h_az_2.^2 + inv_sig_el * h_el_2.^2);
    
    % To save time, we approximate the cost as sum of 1/trace(Inst_FIM) or similar
    % STRICT Liraz approach: trace(inv(Sum FIM)).
    % But we have a cost function per Step. 
    % Let's accumulate elements `a,b,c` over all time steps then invert ONCE.
    
    % NOTE: For the optimizer to work well with `sum` structure, 
    % we usually return a scalar cost.
    % J = trace(inv(Sum_k FIM_k)).
    
    % Since this helper is called for Boat and Target separately:
    FIM_tot = zeros(3,3);
    FIM_tot(1,1) = sum(a);
    FIM_tot(1,2) = sum(b); FIM_tot(2,1) = FIM_tot(1,2);
    FIM_tot(2,2) = sum(c);
    
    % Add epsilon for regularization
    FIM_tot = FIM_tot + 1e-6*eye(3);
    
    J_trace = trace(inv(FIM_tot));
end


%% Dynamics Constraints (Non-Cooperative / Dubins)
function [c, ceq] = nonlcon_3D_Dynamics(X, N_B, N_F, N_tot, Nx, Nu, V_A, g, tf, DBm, DFm, BN, FN)
    
    % Extract Coefficients
    Cx_vec = X(1:Nx*N_tot);
    Cu_vec = X(Nx*N_tot+1:end);
    
    Cx = reshape(Cx_vec, [N_tot, Nx]);
    Cu = reshape(Cu_vec, [N_tot, Nu]);
    
    % Basis Split
    CxB = Cx(1:N_B+1, :); CxF = Cx(N_B+2:end, :);
    CuB = Cu(1:N_B+1, :); CuF = Cu(N_B+2:end, :);
    
    % 1. Reconstruct States & Controls
    Psi_A   = BN * CxB(:,3) + FN * CxF(:,3);
    tan_Phi = BN * CuB(:,1) + FN * CuF(:,1);
    
    % 2. Reconstruct Derivatives
    dX_A   = DBm * CxB(:,1) + DFm * CxF(:,1);
    dY_A   = DBm * CxB(:,2) + DFm * CxF(:,2);
    dPsi_A = DBm * CxB(:,3) + DFm * CxF(:,3);
    
    % 3. Constraints (Dubins)
    % dX = V * cos(Psi) * tf
    % dY = V * sin(Psi) * tf
    % dPsi = (g/V) * tan(Phi) * tf
    
    eq1 = dX_A - tf * V_A * cos(Psi_A);
    eq2 = dY_A - tf * V_A * sin(Psi_A);
    eq3 = dPsi_A - tf * (g/V_A) * tan_Phi;
    
    ceq = [eq1; eq2; eq3];
    c   = []; % Inequality handled in bounds/barrier
end

%% Basis Helper Functions (Identical to Liraz/User's)
function BN = bernsteinMatrix(N, time)
    t = (time - time(1)) / (time(end) - time(1));
    t = t(:);
    BN = zeros(length(t), N + 1);
    for k = 0:N
        log_coeff = gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1);
        term = exp(log_coeff + k .* log(t + 1e-16) + (N - k) .* log(1 - t + 1e-16));
        BN(:, k + 1) = term;
    end
    BN(1, 1) = 1; BN(end, end) = 1; % Endpoint fix
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
    t = time(:); Nt = numel(t); T = time(end);
    if T==0, T=1; end
    w_base = 2*pi*(t/T);
    FN = zeros(Nt, 2*N);
    for k = 1:N
        w = w_base * k;
        FN(:, (k-1)*2 + 1) = cos(w);
        FN(:, (k-1)*2 + 2) = sin(w);
    end
end

function Dm = FourierDiff(N, time)
    t = time(:); Nt = numel(t); T = time(end);
    if T==0, T=1; end
    w_base = 2*pi*(t/T);
    Dm = zeros(Nt, 2*N);
    for k = 1:N
        w = w_base * k;
        rate = k * (2*pi/T);
        Dm(:, (k-1)*2 + 1) = -rate .* sin(w);
        Dm(:, (k-1)*2 + 2) =  rate .* cos(w);
    end
end

function [X_A, Y_A, Psi_A, tan_Phi] = getStates_3D(X, N_B, N_F, N_tot, Nx, Nu, BN, FN)
    Cx = reshape(X(1:Nx*N_tot), [N_tot, Nx]);
    Cu = reshape(X(Nx*N_tot+1:end), [N_tot, Nu]);
    
    CxB = Cx(1:N_B+1, :); CxF = Cx(N_B+2:end, :);
    CuB = Cu(1:N_B+1, :); CuF = Cu(N_B+2:end, :);
    
    X_A   = BN * CxB(:,1) + FN * CxF(:,1);
    Y_A   = BN * CxB(:,2) + FN * CxF(:,2);
    Psi_A = BN * CxB(:,3) + FN * CxF(:,3);
    tan_Phi = BN * CuB(:,1) + FN * CuF(:,1);
end