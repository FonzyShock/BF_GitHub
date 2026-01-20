% Optimal_UAV_Path_RoundTrip_Stable.m
% =========================================================================
% STABLE VERSION: Fixed Singularity via Regularization & Order Reduction
% =========================================================================

clear; clc; close all;

%% ========================================================================
%  1. PROBLEM CONFIGURATION
%  ========================================================================
fprintf('=========================================================================\n');
fprintf('  3D ROUND TRIP - STABLE INITIALIZATION\n');
fprintf('=========================================================================\n\n');

% Mission Parameters
tf          = 100;                      % s, Duration
V_A         = 25;                       % m/s, Speed
z_A         = 100;                      % m, Altitude
g           = 9.81;
control_lim = deg2rad(45);              % Max bank angle
r_NFZ       = 100;                      % No-Fly Zone Radius
r_Com       = 50;                       % Recovery Radius
kappa       = 1e3;                      % Barrier weight

% Locations
x0_T    = [500, 0];                     % Target
x0_B    = [0, 0];                       % Boat (Home)
V_B     = 0;                            

% Basis Configuration (REDUCED ORDER FOR STABILITY)
N_B     = 10;                           % Reduced from 20 -> 10
N_F     = 3;                            % Reduced from 5 -> 3
N_tot   = N_B + 1 + 2*N_F; 
Nx      = 3; % [x, y, psi]
Nu      = 1; % [tan(phi)]

% Generate Basis Matrices
time_opt = linspace(0, 1, N_tot)';
BN  = bernsteinMatrix(N_B, time_opt);
FN  = FourierMatrix(N_F, time_opt);
DBm = bernsteinMatrixDerivative(N_B, time_opt);
DFm = FourierDiff(N_F, time_opt);

Basis   = [BN, FN];
DBasis  = [DBm, DFm];

%% ========================================================================
%  2. STABLE INITIALIZATION (Ridge Regression)
%  ========================================================================
fprintf('Generating Smooth Initial Guess...\n');

% Generate a smooth "Boomerang" shape using pure sinusoids
% This guarantees derivatives exist and are smooth, preventing infinite bank angles.
t_guess = linspace(0, tf, N_tot)';
tau     = linspace(0, pi, N_tot)'; % Half-period for round trip

% Shape: Fly out to 500, fly back to 0. Add width (y) to avoid colliding with NFZ.
dist_T  = x0_T(1); 
x_guess = dist_T * sin(tau);            % 0 -> 500 -> 0
y_guess = -150 * sin(tau).^2;           % Dip away from Y=0 to avoid NFZ at (500,0)

% Calculate Heading (Psi) and Bank (Phi) for this geometric path
dx_g = gradient(x_guess, t_guess);
dy_g = gradient(y_guess, t_guess);
psi_guess = atan2(dy_g, dx_g); 
psi_guess = unwrap(psi_guess);

dpsi_g = gradient(psi_guess, t_guess);
tan_phi_guess = (dpsi_g * V_A) / g;

% Clamp initial guess controls to be valid
max_tan = tan(control_lim);
tan_phi_guess = max(min(tan_phi_guess, max_tan), -max_tan);

% --- RIDGE REGRESSION FIT (The Fix) ---
% Solves (A'A + lambda*I)x = A'b instead of A\b
% This handles singular/ill-conditioned Basis matrices.
lambda = 1e-6; 

Cx_init = zeros(N_tot, Nx);
Cx_init(:,1) = ridge_fit(Basis, x_guess, lambda);
Cx_init(:,2) = ridge_fit(Basis, y_guess, lambda);
Cx_init(:,3) = ridge_fit(Basis, psi_guess, lambda);

Cu_init = zeros(N_tot, Nu);
Cu_init(:,1) = ridge_fit(Basis, tan_phi_guess, lambda);

X0 = [reshape(Cx_init, [Nx*N_tot, 1]); reshape(Cu_init, [Nu*N_tot, 1])];

% Verify Feasibility of Guess
fprintf('Initial Guess Generated. Max X: %.1f, Max Y: %.1f\n', max(x_guess), min(y_guess));

%% ========================================================================
%  3. OPTIMIZATION
%  ========================================================================

% Mission Data
mission.x0_T = x0_T; mission.x0_B = x0_B;
mission.V_B  = V_B;  mission.r_NFZ = r_NFZ;
mission.r_Com = r_Com; mission.z_A = z_A;
mission.kappa = kappa;

% Constraints: Initial Condition (Equality)
Aeq = zeros(Nx, numel(X0)); beq = zeros(Nx, 1);
Basis_0 = [BN(1,:), FN(1,:)]; 
x_start = [x0_B(1), x0_B(2), psi_guess(1)]; 

for i = 1:Nx
    Aeq(i, (i-1)*N_tot+1 : i*N_tot) = Basis_0;
    beq(i) = x_start(i);
end

% Function Handles
fun = @(X) costFunc(X, N_tot, Nx, Nu, mission, tf, Basis);
nonlcon = @(X) constraints(X, N_tot, Nx, Nu, V_A, g, tf, control_lim, Basis, DBasis, mission);

% Options
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
    'MaxFunctionEvaluations', 2e5, 'MaxIterations', 3000, ...
    'ConstraintTolerance', 1e-3, 'OptimalityTolerance', 1e-3, ...
    'StepTolerance', 1e-8);

fprintf('Starting Optimization...\n');
tic;
[xOut, Jout, exitflag] = fmincon(fun, X0, [], [], Aeq, beq, [], [], nonlcon, options);
toc;

%% ========================================================================
%  4. VISUALIZATION
%  ========================================================================
t_plot = linspace(0, 1, 200)';
BN_p = bernsteinMatrix(N_B, t_plot);
FN_p = FourierMatrix(N_F, t_plot);
Basis_p = [BN_p, FN_p];

[X_p, Y_p, Psi_p, tan_Phi_p] = reconstruct(xOut, N_tot, Nx, Nu, Basis_p);

figure('Color','w','Position',[100 100 1000 500]);

% Map
subplot(1,2,1); hold on; grid on; axis equal;
plot(X_p, Y_p, 'b-', 'LineWidth', 2, 'DisplayName', 'UAV Path');
plot(x0_B(1), x0_B(2), 'gs', 'MarkerFaceColor','g', 'DisplayName', 'Boat');
plot(x0_T(1), x0_T(2), 'rd', 'MarkerFaceColor','r', 'DisplayName', 'Target');
viscircles(x0_T, r_NFZ, 'Color', 'r', 'LineStyle', '--');
viscircles(x0_B, r_Com, 'Color', 'g', 'LineStyle', '-.');
legend('Location','best'); xlabel('X'); ylabel('Y'); title('Optimal Path');

% Controls
subplot(1,2,2); hold on; grid on;
phi_deg = rad2deg(atan(tan_Phi_p));
plot(t_plot*tf, phi_deg, 'LineWidth', 2);
yline(rad2deg(control_lim), 'r--'); yline(-rad2deg(control_lim), 'r--');
xlabel('Time [s]'); ylabel('Bank Angle [deg]'); title('Control Inputs');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

% Robust Fitting Function
function c = ridge_fit(A, b, lambda)
    % (A'A + lambda*I) * c = A'b
    c = (A'*A + lambda*eye(size(A,2))) \ (A'*b);
end

function J = costFunc(X, N_tot, Nx, Nu, data, tf, Basis)
    [X_A, Y_A, ~, ~] = reconstruct(X, N_tot, Nx, Nu, Basis);
    
    % Info Gain Proxy: Minimize average squared distance to target
    dist_sq = (X_A - data.x0_T(1)).^2 + (Y_A - data.x0_T(2)).^2;
    dist = sqrt(dist_sq);
    
    % Soft Barrier for NFZ
    violation = max(0, data.r_NFZ - dist);
    J_barrier = data.kappa * sum(violation.^2);
    
    J = sum(dist_sq) + J_barrier; 
end

function [c, ceq] = constraints(X, N_tot, Nx, Nu, V_A, g, tf, lim, Basis, DBasis, data)
    [X_A, Y_A, Psi_A, tan_Phi] = reconstruct(X, N_tot, Nx, Nu, Basis);
    
    % Derivatives via Matrix Multiplication
    Cx = reshape(X(1:Nx*N_tot), [N_tot, Nx]);
    dX_A   = DBasis * Cx(:,1);
    dY_A   = DBasis * Cx(:,2);
    dPsi_A = DBasis * Cx(:,3);
    
    % Dynamics (Dubins)
    eq1 = dX_A - tf * V_A * cos(Psi_A);
    eq2 = dY_A - tf * V_A * sin(Psi_A);
    eq3 = dPsi_A - tf * (g/V_A) * tan_Phi;
    ceq = [eq1; eq2; eq3];
    
    % Inequality
    % 1. Control Limits
    c_ctrl = [tan_Phi - tan(lim); -tan(lim) - tan_Phi];
    
    % 2. Terminal Constraint (Return to Boat)
    final_dist_sq = (X_A(end) - data.x0_B(1))^2 + (Y_A(end) - data.x0_B(2))^2;
    c_term = final_dist_sq - data.r_Com^2;
    
    c = [c_ctrl; c_term];
end

function [X_A, Y_A, Psi_A, tan_Phi] = reconstruct(X, N_tot, Nx, Nu, Basis)
    Cx = reshape(X(1:Nx*N_tot), [N_tot, Nx]);
    Cu = reshape(X(Nx*N_tot+1:end), [N_tot, Nu]);
    X_A   = Basis * Cx(:,1);
    Y_A   = Basis * Cx(:,2);
    Psi_A = Basis * Cx(:,3);
    tan_Phi = Basis * Cu(:,1);
end

% Standard Basis
function BN = bernsteinMatrix(N, t)
    t=t(:); BN=zeros(length(t),N+1);
    for k=0:N, BN(:,k+1) = nchoosek(N,k).*(t.^k).*((1-t).^(N-k)); end
end
function BN_dot = bernsteinMatrixDerivative(N, t)
    t=t(:); BN_m1 = bernsteinMatrix(N-1, t);
    BN_dot = N * ([zeros(length(t),1), BN_m1] - [BN_m1, zeros(length(t),1)]);
end
function FN = FourierMatrix(N, t)
    t=t(:); FN=zeros(length(t),2*N);
    for k=1:N, w=2*pi*k*t; FN(:,2*k-1)=cos(w); FN(:,2*k)=sin(w); end
end
function DF = FourierDiff(N, t)
    t=t(:); DF=zeros(length(t),2*N);
    for k=1:N, w=2*pi*k*t; r=2*pi*k; DF(:,2*k-1)=-r*sin(w); DF(:,2*k)=r*cos(w); end
end