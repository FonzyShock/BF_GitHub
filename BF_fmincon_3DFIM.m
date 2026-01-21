% moving_target_observer_3D_FIM.m
% =========================================================================
% Optimal Observer Trajectory Generation for Moving Target Localization
% Using Mixed Bernstein-Fourier Approximants
%
% Updates:
%   - Cost function updated to use Accumulated 3D Fisher Information Matrix
%
% Author: Liraz Mudrik (Modified for 3D FIM)
% =========================================================================
clear; clc; close all;

%% ========================================================================
%  COMMON PROBLEM PARAMETERS
%  ========================================================================
fprintf('=========================================================================\n');
fprintf('  MOVING TARGET OBSERVER TRAJECTORY OPTIMIZATION (3D FIM)\n');
fprintf('=========================================================================\n\n');

% UAV parameters
V_A   = 30;                      % m/s, UAV constant velocity
z_A   = 100;                     % m, UAV altitude
psi0  = 0;                       % rad, UAV initial heading
x0_A  = [0, 0];                  % m, UAV initial position
u0    = [cos(psi0), sin(psi0)];  % initial control (heading direction)

% Target parameters (MOVING TARGET)
x_T0  = [400, 0];                % m, target initial position
v_T   = [5, 0];                  % m/s, target velocity (known)

% Measurement parameters
Sigma = diag([3e-3, 3e-3].^2);   % Measurement noise covariance (rad^2)
freq  = 10;                      % Hz, measurement frequency

% Constraint parameters
r_NFZ = 100;                     % m, no-fly zone radius
kappa = 1e-1;                    % Barrier parameter (tunable)

% Time parameters
tf = 90;                         % s, scenario duration

% Approximation orders
N_B = 60;                        % Bernstein polynomial order
N_F = 24;                        % Fourier series order

% State and control dimensions
Nx = 2;                          % number of states (x, y)
Nu = 2;                          % number of controls (cos(psi), sin(psi))

%% ========================================================================
%  PART 1: MIXED BERNSTEIN-FOURIER WITH MOVING TARGET
%  ========================================================================
fprintf('------------------------------------------------------------------\n');
fprintf('  PART 1: Mixed Bernstein-Fourier Optimization (Moving Target)\n');
fprintf('------------------------------------------------------------------\n');

% Total number of coefficients
N_tot = N_B + 1 + 2*N_F + 1;

% Time nodes
Nt_BF = N_B + 1;
time_BF = linspace(0, 1, Nt_BF);

% Construct basis matrices
BN_BF  = bernsteinMatrix(N_B, time_BF);
FN_BF  = FourierMatrix(N_F, time_BF);
DBm_BF = bernsteinMatrixDerivative(N_B, time_BF);
DFm_BF = FourierDiff(N_F, time_BF);

fprintf('Basis matrix condition number: %.2e\n', cond([BN_BF, FN_BF]));

% Initial guess: fly toward target, then orbit
Cx_BF = zeros(N_tot, Nx);
Cu_BF = zeros(N_tot, Nu);

for k = 1:Nt_BF
    t_k = time_BF(k) * tf;
    x_T_k = x_T0(1) + v_T(1) * t_k;
    y_T_k = x_T0(2) + v_T(2) * t_k;
    
    t_approach = norm(x_T0 + [-r_NFZ, 0]) / V_A;
    
    if t_k <= t_approach
        Cx_BF(k, :) = [V_A * t_k, 0];
        Cu_BF(k, :) = [1, 0];
    else
        angle = pi - V_A / r_NFZ * (t_k - t_approach);
        Cx_BF(k, :) = [x_T_k, y_T_k] + r_NFZ * [cos(angle), sin(angle)];
        psi_k = pi/2 - V_A / r_NFZ * (t_k - t_approach);
        Cu_BF(k, :) = [cos(psi_k), sin(psi_k)];
    end
end

% Reshape for optimizer
X0_BF = [reshape(Cx_BF, [Nx * N_tot, 1]); reshape(Cu_BF, [Nu * N_tot, 1])];

% Equality constraints (initial conditions)
Aeq_BF = blkdiag([BN_BF(1,:), FN_BF(1,:)], [BN_BF(1,:), FN_BF(1,:)], ...
                 [BN_BF(1,:), FN_BF(1,:)], [BN_BF(1,:), FN_BF(1,:)]);
beq_BF = [x0_A'; u0'];

% Bounds
lb_BF = -inf(size(X0_BF));
ub_BF =  inf(size(X0_BF));

% Optimizer options
options_BF = optimoptions(@fmincon, ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1e6, ...
    'OptimalityTolerance', 1e-4,...
    'ConstraintTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'MaxIterations', 5e3, ...
    'Display', 'iter');

% Define cost and constraint functions
costFunc_BF = @(X) costFunc_MixedBF_MovingTarget_3DFIM(X, N_B, N_F, N_tot, Nx, Nu, ...
    x_T0, v_T, z_A, Sigma, tf, r_NFZ, kappa, BN_BF, FN_BF, time_BF);

nonlcon_BF = @(X) nonlcon_MixedBF_MovingTarget(X, N_B, N_F, N_tot, Nx, Nu, ...
    V_A, tf, DBm_BF, DFm_BF, BN_BF, FN_BF);

% Run optimization (cold start)
fprintf('\nStarting Mixed BF cold start optimization...\n');
tic;
[xOut_BF_cold, Jout_BF_cold, exitflag_BF_cold, output_BF_cold, lambda_BF_cold] = ...
    fmincon(costFunc_BF, X0_BF, [], [], Aeq_BF, beq_BF, lb_BF, ub_BF, nonlcon_BF, options_BF);
time_BF_cold = toc;

if exitflag_BF_cold < 1 && isfield(output_BF_cold, 'bestfeasible') && ~isempty(output_BF_cold.bestfeasible)
    xOut_BF_cold = output_BF_cold.bestfeasible.x;
    Jout_BF_cold = output_BF_cold.bestfeasible.fval;
    fprintf('Using best feasible solution for cold start.\n');
end

fprintf('Mixed BF cold start completed in %.2f seconds\n', time_BF_cold);
fprintf('  Cost (CRLB + barrier): %.5e\n', Jout_BF_cold);
fprintf('  Exit flag: %d\n', exitflag_BF_cold);
fprintf('  Iterations: %d\n', output_BF_cold.iterations);

% Store Mixed BF results
results.mixedBF.xOut = xOut_BF_cold;
results.mixedBF.lambda = lambda_BF_cold;
results.mixedBF.time_cold = time_BF_cold;
results.mixedBF.iterations_cold = output_BF_cold.iterations;
results.mixedBF.exitflag = exitflag_BF_cold;
results.mixedBF.numVars = length(X0_BF);

%% ========================================================================
%  PART 2: BERNSTEIN-ONLY WITH MOVING TARGET (BASELINE)
%  ========================================================================
fprintf('\n------------------------------------------------------------------\n');
fprintf('  PART 2: Bernstein-Only Optimization (Moving Target, Baseline)\n');
fprintf('------------------------------------------------------------------\n');

N_bern = 500;
r_NFZ_bern = 10; % Relaxed constraint for baseline

% Time nodes
Nt_bern = N_bern + 1;
time_bern = linspace(0, 1, Nt_bern);

% Basis matrices
BN_bern = bernsteinMatrix(N_bern, time_bern);
Dm_bern = Diff_elev(N_bern, 1);

% Initial guess
Cx_bern = zeros(N_bern + 1, Nx);
Cu_bern = zeros(N_bern + 1, Nu);
for k = 1:Nt_bern
    t_k = time_bern(k) * tf;
    x_T_k = x_T0(1) + v_T(1) * t_k;
    y_T_k = x_T0(2) + v_T(2) * t_k;
    
    t_approach = norm(x_T0 + [-r_NFZ_bern, 0]) / V_A;
    
    if t_k <= t_approach
        Cx_bern(k, :) = [V_A * t_k, 0];
        Cu_bern(k, :) = [1, 0];
    else
        angle = pi - V_A / r_NFZ_bern * (t_k - t_approach);
        Cx_bern(k, :) = [x_T_k, y_T_k] + r_NFZ_bern * [cos(angle), sin(angle)];
        psi_k = pi/2 - V_A / r_NFZ_bern * (t_k - t_approach);
        Cu_bern(k, :) = [cos(psi_k), sin(psi_k)];
    end
end
X0_bern = [reshape(Cx_bern, [Nx * (N_bern + 1), 1]); ...
           reshape(Cu_bern, [Nu * (N_bern + 1), 1])];

% Equality constraints
Aeq_bern = zeros(Nx + Nu, (N_bern + 1) * (Nx + Nu));
Aeq_bern(1, 1) = 1;
Aeq_bern(2, N_bern + 2) = 1;
Aeq_bern(3, (N_bern + 1) * Nx + 1) = 1;
Aeq_bern(4, (N_bern + 1) * Nx + N_bern + 2) = 1;
beq_bern = [x0_A'; u0'];

% Bounds
lb_bern = -inf(size(X0_bern));
ub_bern =  inf(size(X0_bern));

% Optimizer options
options_bern = optimoptions(@fmincon, ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1e6, ...
    'ConstraintTolerance', 1e-6, ...
    'StepTolerance', 1e-8, ...
    'MaxIterations', 5e3, ...
    'Display', 'iter');

% Define cost and constraint functions
costFunc_bern = @(X) costFunc_Bernstein_MovingTarget_3DFIM(X, N_bern, Nx, Nu, ...
    x_T0, v_T, z_A, Sigma, tf, r_NFZ, kappa, time_bern);
nonlcon_bern = @(X) nonlcon_Bernstein_MovingTarget(X, N_bern, Nx, Nu, ...
    V_A, tf, r_NFZ_bern, x_T0, v_T, Dm_bern, time_bern);

% Run optimization
fprintf('\nStarting Bernstein-only optimization...\n');
tic;
[xOut_bern, Jout_bern, exitflag_bern, output_bern, lambda_bern] = ...
    fmincon(costFunc_bern, X0_bern, [], [], Aeq_bern, beq_bern, lb_bern, ub_bern, nonlcon_bern, options_bern);
time_bern_elapsed = toc;

fprintf('Bernstein-only completed in %.2f seconds\n', time_bern_elapsed);
fprintf('  Cost (CRLB): %.5e\n', Jout_bern);

% Store Bernstein results
results.bernstein.xOut = xOut_bern;
results.bernstein.lambda = lambda_bern;
results.bernstein.time = time_bern_elapsed;
results.bernstein.iterations = output_bern.iterations;
results.bernstein.exitflag = exitflag_bern;
results.bernstein.numVars = length(X0_bern);

%% ========================================================================
%  PART 3: EVALUATE AND COMPARE RESULTS
%  ========================================================================
fprintf('\n------------------------------------------------------------------\n');
fprintf('  PART 3: Evaluate and Compare Results\n');
fprintf('------------------------------------------------------------------\n');

% Fine evaluation grid (10 Hz)
Nt_eval = tf * freq + 1;
time_eval = linspace(0, 1, Nt_eval);
time_eval_sec = time_eval * tf;

% Evaluation basis matrices
BN_eval = bernsteinMatrix(N_B, time_eval);
FN_eval = FourierMatrix(N_F, time_eval);
BN_eval_bern = bernsteinMatrix(N_bern, time_eval);

% Extract Mixed BF trajectories
[X_A_BF, Y_A_BF, cos_U_BF, sin_U_BF] = ...
    getFullStates_MixedBF(results.mixedBF.xOut, N_B, N_F, N_tot, Nx, Nu, BN_eval, FN_eval);

% Extract Bernstein trajectories
[X_A_bern, Y_A_bern] = getStates_Bernstein(results.bernstein.xOut, N_bern, Nx, BN_eval_bern);

% Compute target trajectory
X_T_traj = x_T0(1) + v_T(1) * time_eval_sec';
Y_T_traj = x_T0(2) + v_T(2) * time_eval_sec';

% Compute distances to target
dist_BF = sqrt((X_A_BF - X_T_traj).^2 + (Y_A_BF - Y_T_traj).^2);
dist_bern = sqrt((X_A_bern - X_T_traj).^2 + (Y_A_bern - Y_T_traj).^2);

% Compute CRLB trace for both methods (Using 3D FIM accumulation)
[CRLB_BF, CRLB_BF_cumulative] = computeCRLB_MovingTarget_3DFIM(X_A_BF, Y_A_BF, x_T0, v_T, z_A, Sigma, tf, time_eval);
[CRLB_bern, CRLB_bern_cumulative] = computeCRLB_MovingTarget_3DFIM(X_A_bern, Y_A_bern, x_T0, v_T, z_A, Sigma, tf, time_eval);

% Store evaluation results
results.eval.time = time_eval_sec;
results.eval.X_A_BF = X_A_BF;
results.eval.Y_A_BF = Y_A_BF;
results.eval.X_A_bern = X_A_bern;
results.eval.Y_A_bern = Y_A_bern;
results.eval.X_T = X_T_traj;
results.eval.Y_T = Y_T_traj;
results.eval.dist_BF = dist_BF;
results.eval.dist_bern = dist_bern;
results.eval.CRLB_BF = CRLB_BF;
results.eval.CRLB_bern = CRLB_bern;
results.eval.CRLB_BF_cumulative = CRLB_BF_cumulative;
results.eval.CRLB_bern_cumulative = CRLB_bern_cumulative;

fprintf('\n--- Performance Comparison ---\n');
fprintf('%-30s %15s %15s\n', 'Metric', 'Mixed BF', 'Bernstein');
fprintf('%-30s %15s %15s\n', repmat('-',1,30), repmat('-',1,15), repmat('-',1,15));
fprintf('%-30s %15.5e %15.5e\n', 'Final CRLB trace', CRLB_BF, CRLB_bern);
fprintf('%-30s %15.2f %15.2f\n', 'Computation time (s)', results.mixedBF.time_cold, results.bernstein.time);
fprintf('%-30s %15d %15d\n', 'Decision variables', results.mixedBF.numVars, results.bernstein.numVars);

CRLB_improvement = (CRLB_bern - CRLB_BF) / CRLB_bern * 100;
var_reduction = (results.bernstein.numVars - results.mixedBF.numVars) / results.bernstein.numVars * 100;

fprintf('\n--- Improvement Metrics ---\n');
fprintf('CRLB improvement (Mixed BF vs Bernstein): %.2f%%\n', CRLB_improvement);
fprintf('Variable reduction: %.2f%%\n', var_reduction);

%% ========================================================================
%  PART 4: COVECTOR MAPPING VALIDATION
%  ========================================================================
fprintf('\n------------------------------------------------------------------\n');
fprintf('  PART 4: Covector Mapping Validation\n');
fprintf('------------------------------------------------------------------\n');

Nt_opt = Nt_BF;
if isfield(results.mixedBF.lambda, 'eqnonlin') && ~isempty(results.mixedBF.lambda.eqnonlin)
    lambda_eq = results.mixedBF.lambda.eqnonlin;
    lambda_x_NLP = lambda_eq(1:Nt_opt);
    lambda_y_NLP = lambda_eq(Nt_opt+1:2*Nt_opt);
    
    w = 1 / (Nt_opt + 1);
    lambda_x_scaled = lambda_x_NLP / w;
    lambda_y_scaled = lambda_y_NLP / w;
    
    results.covector.lambda_x_NLP = lambda_x_scaled;
    results.covector.lambda_y_NLP = lambda_y_scaled;
else
    results.covector.lambda_x_NLP = [];
    results.covector.lambda_y_NLP = [];
end

[H_traj, H_components] = computeHamiltonian_MovingTarget_3DFIM(...
    X_A_BF, Y_A_BF, cos_U_BF, sin_U_BF, ...
    x_T0, v_T, z_A, Sigma, tf, r_NFZ, kappa, V_A, ...
    results.covector.lambda_x_NLP, results.covector.lambda_y_NLP, ...
    time_eval, BN_eval, FN_eval, N_B, N_F);

results.covector.H = H_traj;
results.covector.H_components = H_components;

fprintf('\nHamiltonian analysis (low-order, N_B=%d, N_F=%d):\n', N_B, N_F);
fprintf('  Mean: %.6f\n', mean(H_traj));
fprintf('  Std:  %.6f\n', std(H_traj));

%% ========================================================================
%  PART 6: GENERATE FIGURES
%  ========================================================================
fprintf('\n------------------------------------------------------------------\n');
fprintf('  PART 6: Generate Figures\n');
fprintf('------------------------------------------------------------------\n');

% Figure 5: Top-down trajectories
figure('Name', 'Fig5_TopDown_Trajectories', 'Position', [100, 100, 800, 600]);
hold on; grid on; axis equal;
plot(X_A_BF, Y_A_BF, 'b-', 'LineWidth', 2, 'DisplayName', 'Mixed BF');
plot(X_A_bern, Y_A_bern, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Bernstein-only');
plot(X_T_traj, Y_T_traj, 'r--', 'LineWidth', 2, 'DisplayName', 'Target');
plot(x_T0(1), x_T0(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Target start');
plot(X_T_traj(end), Y_T_traj(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Target end');

theta_circle = linspace(0, 2*pi, 100);
plot(x_T0(1) + r_NFZ * cos(theta_circle), x_T0(2) + r_NFZ * sin(theta_circle), 'r:', 'LineWidth', 1.5);
plot(X_T_traj(end) + r_NFZ * cos(theta_circle), Y_T_traj(end) + r_NFZ * sin(theta_circle), 'm:', 'LineWidth', 1.5, 'DisplayName', 'NFZ (End)');
plot(x0_A(1), x0_A(2), 'b^', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'UAV start');
xlabel('x, m'); ylabel('y, m');
legend('Location', 'best');
title(sprintf('Optimal UAV trajectories (3D FIM Trace: BF=%.2e, Bern=%.2e)', CRLB_BF, CRLB_bern));

% Figure 7: Constraint and information
figure('Name', 'Fig7_Constraint_Information', 'Position', [100, 100, 800, 500]);
subplot(2, 1, 1);
hold on; grid on;
plot(time_eval_sec, dist_BF, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Mixed BF');
plot(time_eval_sec, dist_bern, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Bernstein-only');
yline(r_NFZ, 'r--', 'NFZ');
xlabel('Time, s'); ylabel('Distance to target, m');
legend('Location', 'best');
title('UAV-Target distance');

subplot(2, 1, 2);
hold on; grid on;
% Plot cumulative trace of inverse FIM
plot(time_eval_sec, CRLB_BF_cumulative, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Mixed BF');
plot(time_eval_sec, CRLB_bern_cumulative, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Bernstein-only');
xlabel('Time, s'); ylabel('Trace(FIM^{-1})');
legend('Location', 'best');
title('Cumulative A-Optimality Cost');

% Save results
save('moving_target_results.mat', 'results');

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

%% Cost function for Mixed BF (3D FIM Accumulated)
function J = costFunc_MixedBF_MovingTarget_3DFIM(X, N_B, N_F, N_tot, Nx, Nu, ...
    x_T0, v_T, z_A, Sigma, tf, r_NFZ, kappa, BN, FN, time)
    %#ok<INUSD> 
    Cx = reshape(X(1:Nx*N_tot), [N_tot, Nx]);
    
    % States on nodes
    CxB = Cx(1:N_B+1, :);
    CxF = Cx(N_B+2:end, :);
    X_A = BN * CxB(:,1) + FN * CxF(:,1);
    Y_A = BN * CxB(:,2) + FN * CxF(:,2);
    
    % Target trajectory on nodes
    tsec = (time(:) * tf);
    X_T  = x_T0(1) + v_T(1) * tsec;
    Y_T  = x_T0(2) + v_T(2) * tsec;
    Z_T  = zeros(size(X_T)); 
    
    % 3D Relative Geometry (r = Target - Cam)
    rx = X_T - X_A;
    ry = Y_T - Y_A;
    rz = Z_T - z_A; % This will be negative altitude
    
    eps_guard = 1e-6;
    rxy2 = rx.^2 + ry.^2 + eps_guard;
    rxy  = sqrt(rxy2);
    r2   = rxy2 + rz.^2 + eps_guard;
    
    % Jacobian Components
    % H_az = [ -ry ./ rxy2;  rx ./ rxy2;  0 ]
    h_az_x = -ry ./ rxy2;
    h_az_y =  rx ./ rxy2;
    h_az_z =  zeros(size(rx));
    
    % H_el = [ (rx.*rz)./(r2.*rxy); (ry.*rz)./(r2.*rxy); -rxy./r2 ]
    h_el_x = (rx .* rz) ./ (r2 .* rxy);
    h_el_y = (ry .* rz) ./ (r2 .* rxy);
    h_el_z = -rxy ./ r2;
    
    % Weights
    w_az = 1 / Sigma(1,1);
    w_el = 1 / Sigma(2,2);
    
    % Accumulate FIM (Vectorized Summation)
    F11 = w_az * sum(h_az_x.^2) + w_el * sum(h_el_x.^2);
    F22 = w_az * sum(h_az_y.^2) + w_el * sum(h_el_y.^2);
    F33 = w_az * sum(h_az_z.^2) + w_el * sum(h_el_z.^2);
    F12 = w_az * sum(h_az_x .* h_az_y) + w_el * sum(h_el_x .* h_el_y);
    F13 = w_az * sum(h_az_x .* h_az_z) + w_el * sum(h_el_x .* h_el_z);
    F23 = w_az * sum(h_az_y .* h_az_z) + w_el * sum(h_el_y .* h_el_z);
    
    FIM = [F11, F12, F13;
           F12, F22, F23;
           F13, F23, F33];
       
    % A-Optimality (Trace of Inverse)
    % Add regularization for numerical stability of inverse
    CRLB_trace = trace(inv(FIM + 1e-9*eye(3)));
    
    % Barrier Function (NFZ is defined in 2D plane usually)
    % r_xy distance used for NFZ
    barrier_term = kappa * norm(max(0, r_NFZ - rxy))^2;
    
    J = CRLB_trace + barrier_term;
end

%% Cost function for Bernstein (3D FIM Accumulated)
function J = costFunc_Bernstein_MovingTarget_3DFIM(X, N, Nx, Nu, x_T0, v_T, z_A, Sigma, tf, r_NFZ, kappa, time)
    
    Cx = reshape(X(1:Nx*(N+1)), [N+1, Nx]);
    X_A = Cx(:, 1);
    Y_A = Cx(:, 2);
    
    tsec = time(:) * tf;
    X_T = x_T0(1) + v_T(1) * tsec;
    Y_T = x_T0(2) + v_T(2) * tsec;
    Z_T = zeros(size(X_T));
    
    rx = X_T - X_A;
    ry = Y_T - Y_A;
    rz = Z_T - z_A;
    
    eps_guard = 1e-6;
    rxy2 = rx.^2 + ry.^2 + eps_guard;
    rxy  = sqrt(rxy2);
    r2   = rxy2 + rz.^2 + eps_guard;
    
    h_az_x = -ry ./ rxy2;
    h_az_y =  rx ./ rxy2;
    h_az_z =  zeros(size(rx));
    
    h_el_x = (rx .* rz) ./ (r2 .* rxy);
    h_el_y = (ry .* rz) ./ (r2 .* rxy);
    h_el_z = -rxy ./ r2;
    
    w_az = 1 / Sigma(1,1);
    w_el = 1 / Sigma(2,2);
    
    F11 = w_az * sum(h_az_x.^2) + w_el * sum(h_el_x.^2);
    F22 = w_az * sum(h_az_y.^2) + w_el * sum(h_el_y.^2);
    F33 = w_az * sum(h_az_z.^2) + w_el * sum(h_el_z.^2);
    F12 = w_az * sum(h_az_x .* h_az_y) + w_el * sum(h_el_x .* h_el_y);
    F13 = w_az * sum(h_az_x .* h_az_z) + w_el * sum(h_el_x .* h_el_z);
    F23 = w_az * sum(h_az_y .* h_az_z) + w_el * sum(h_el_y .* h_el_z);
    
    FIM = [F11, F12, F13; F12, F22, F23; F13, F23, F33];
    CRLB_trace = trace(inv(FIM + 1e-9*eye(3)));
    
    barrier_term = kappa * norm(max(0, r_NFZ - rxy))^2;
    J = CRLB_trace + barrier_term;
end

%% Compute CRLB History (3D FIM Accumulated)
function [CRLB_final, CRLB_cumulative] = computeCRLB_MovingTarget_3DFIM(X_A, Y_A, x_T0, v_T, z_A, Sigma, tf, time)
    
    Nt = length(time);
    tsec = time * tf;
    
    X_T = x_T0(1) + v_T(1) * tsec(:);
    Y_T = x_T0(2) + v_T(2) * tsec(:);
    Z_T = zeros(Nt, 1);
    
    rx = X_T - X_A;
    ry = Y_T - Y_A;
    rz = Z_T - z_A;
    
    eps_guard = 1e-6;
    rxy2 = rx.^2 + ry.^2 + eps_guard;
    rxy  = sqrt(rxy2);
    r2   = rxy2 + rz.^2 + eps_guard;
    
    h_az_x = -ry ./ rxy2;
    h_az_y =  rx ./ rxy2;
    h_az_z =  zeros(size(rx));
    
    h_el_x = (rx .* rz) ./ (r2 .* rxy);
    h_el_y = (ry .* rz) ./ (r2 .* rxy);
    h_el_z = -rxy ./ r2;
    
    w_az = 1 / Sigma(1,1);
    w_el = 1 / Sigma(2,2);
    
    CRLB_cumulative = zeros(Nt, 1);
    
    % Initialize running sums of FIM elements
    f11_sum = 0; f22_sum = 0; f33_sum = 0;
    f12_sum = 0; f13_sum = 0; f23_sum = 0;
    
    for k = 1:Nt
        % Add contribution of current time step k
        f11_sum = f11_sum + w_az * h_az_x(k)^2 + w_el * h_el_x(k)^2;
        f22_sum = f22_sum + w_az * h_az_y(k)^2 + w_el * h_el_y(k)^2;
        f33_sum = f33_sum + w_az * h_az_z(k)^2 + w_el * h_el_z(k)^2;
        
        f12_sum = f12_sum + w_az * h_az_x(k)*h_az_y(k) + w_el * h_el_x(k)*h_el_y(k);
        f13_sum = f13_sum + w_az * h_az_x(k)*h_az_z(k) + w_el * h_el_x(k)*h_el_z(k);
        f23_sum = f23_sum + w_az * h_az_y(k)*h_az_z(k) + w_el * h_el_y(k)*h_el_z(k);
        
        FIM_k = [f11_sum, f12_sum, f13_sum;
                 f12_sum, f22_sum, f23_sum;
                 f13_sum, f23_sum, f33_sum];
             
        if k < 5 || rcond(FIM_k) < 1e-12
            CRLB_cumulative(k) = NaN; % Not observable yet
        else
            CRLB_cumulative(k) = trace(inv(FIM_k));
        end
    end
    
    CRLB_final = CRLB_cumulative(end);
end

%% Compute Hamiltonian (3D FIM Incremental contribution)
function [H, H_components] = computeHamiltonian_MovingTarget_3DFIM(...
    X_A, Y_A, cos_U, sin_U, x_T0, v_T, z_A, Sigma, tf, r_NFZ, kappa, V_A, ...
    lambda_x, lambda_y, time, BN, FN, N_B, N_F)
    
    Nt = length(X_A);
    tsec = time * tf;
    
    X_T = x_T0(1) + v_T(1) * tsec(:);
    Y_T = x_T0(2) + v_T(2) * tsec(:);
    
    rx = X_T - X_A;
    ry = Y_T - Y_A;
    rz = -z_A;
    
    eps_guard = 1e-6;
    rxy2 = rx.^2 + ry.^2 + eps_guard;
    rxy  = sqrt(rxy2);
    r2   = rxy2 + rz.^2 + eps_guard;
    
    % Hamiltonian for an integral cost J = Integral(L dt)
    % Here our cost is J = Trace(Inverse(Sum(FIM_k)))
    % This is a terminal cost on an accumulated state, not a standard integral cost.
    % However, to check costates, we treat the incremental FIM magnitude as a proxy for "L"
    
    L_FIM_Proxy = zeros(Nt, 1);
    w_az = 1/Sigma(1,1); w_el = 1/Sigma(2,2);
    
    for k = 1:Nt
         % Magnitude of information gain at this step
         ha_x = -ry(k)/rxy2(k); ha_y = rx(k)/rxy2(k);
         he_x = rx(k)*rz/ (r2(k)*rxy(k)); he_y = ry(k)*rz/(r2(k)*rxy(k)); he_z = -rxy(k)/r2(k);
         
         % Trace of the instantaneous FIM (Scalar representing info rate)
         tr_FIM_inst = w_az*(ha_x^2 + ha_y^2) + w_el*(he_x^2 + he_y^2 + he_z^2);
         L_FIM_Proxy(k) = tr_FIM_inst;
    end
    
    % Barrier
    L_barrier = kappa * norm(max(0,r_NFZ - rxy))^2;
    
    % Dynamics
    Nt_costate = length(lambda_x);
    time_costate = linspace(0, 1, Nt_costate);
    lambda_x_interp = interp1(time_costate, lambda_x, time, 'linear', 'extrap');
    lambda_y_interp = interp1(time_costate, lambda_y, time, 'linear', 'extrap');
    
    H_dynamics = lambda_x_interp(:) .* (V_A * cos_U) + lambda_y_interp(:) .* (V_A * sin_U);
    
    H = L_FIM_Proxy + L_barrier + H_dynamics / tf;
    
    H_components.L_FIM = L_FIM_Proxy;
    H_components.L_barrier = L_barrier;
    H_components.H_dynamics = H_dynamics;
end

%% Nonlinear constraints for Mixed BF with moving target (dynamics only)
function [c, ceq] = nonlcon_MixedBF_MovingTarget(X, N_B, N_F, N_tot, Nx, Nu, ...
    V_A, tf, DBm, DFm, BN, FN)
    
    Cx = reshape(X(1:Nx*N_tot), [N_tot, Nx]);
    Cu = reshape(X(Nx*N_tot+1:end), [N_tot, Nu]);
    
    % Control inputs on nodes
    cos_U = BN * Cu(1:N_B+1, 1) + FN * Cu(N_B+2:end, 1);
    sin_U = BN * Cu(1:N_B+1, 2) + FN * Cu(N_B+2:end, 2);
    
    % State derivatives on nodes
    dX_A = DBm * Cx(1:N_B+1, 1) + DFm * Cx(N_B+2:end, 1);
    dY_A = DBm * Cx(1:N_B+1, 2) + DFm * Cx(N_B+2:end, 2);
    
    % Dynamics constraints
    ceq = [dX_A - tf * V_A * cos_U; 
           dY_A - tf * V_A * sin_U; 
           cos_U.^2 + sin_U.^2 - 1];
    
    c = [];
end

%% Nonlinear constraints for Bernstein
function [c, ceq] = nonlcon_Bernstein_MovingTarget(X, N, Nx, Nu, V_A, tf, r_NFZ, x_T0, v_T, Dm, time)
    Cx = reshape(X(1:Nx*(N+1)), [N+1, Nx]);
    Cu = reshape(X(Nx*(N+1)+1:end), [N+1, Nu]);
    
    cos_U = Cu(:, 1);
    sin_U = Cu(:, 2);
    
    dX_A = Dm' * Cx(:, 1);
    dY_A = Dm' * Cx(:, 2);
    
    ceq = [dX_A - tf * V_A * cos_U; 
           dY_A - tf * V_A * sin_U; 
           cos_U.^2 + sin_U.^2 - 1];
    c = [];
end

%% Get full states for Mixed BF
function [X_A, Y_A, cos_U, sin_U] = getFullStates_MixedBF(X, N_B, N_F, N_tot, Nx, Nu, BN, FN)
    Cx = reshape(X(1:Nx*N_tot), [N_tot, Nx]);
    Cu = reshape(X(Nx*N_tot+1:end), [N_tot, Nu]);
    
    X_A = BN * Cx(1:N_B+1, 1) + FN * Cx(N_B+2:end, 1);
    Y_A = BN * Cx(1:N_B+1, 2) + FN * Cx(N_B+2:end, 2);
    cos_U = BN * Cu(1:N_B+1, 1) + FN * Cu(N_B+2:end, 1);
    sin_U = BN * Cu(1:N_B+1, 2) + FN * Cu(N_B+2:end, 2);
end

%% Get states for Bernstein
function [X_A, Y_A] = getStates_Bernstein(X, N, Nx, BN)
    Cx = reshape(X(1:Nx*(N+1)), [N+1, Nx]);
    X_A = BN * Cx(:, 1);
    Y_A = BN * Cx(:, 2);
end

%% Bernstein matrix
function BN = bernsteinMatrix(N, time)
    t = (time - time(1)) / (time(end) - time(1));
    len_t = length(t);
    BN = zeros(len_t, N + 1);
    
    idx0 = (t == 0);
    idx1 = (t == 1);
    idxInterior = ~(idx0 | idx1);
    
    logt = log(t(idxInterior));
    log1mt = log(1 - t(idxInterior));
    
    for k = 0:N
        log_coeff = gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1);
        BN(idxInterior, k + 1) = exp(log_coeff + k .* logt + (N - k) .* log1mt);
    end
    
    BN(idx0, 1) = 1;
    BN(idx1, N + 1) = 1;
end
%% Bernstein derivative matrix
function BN_dot = bernsteinMatrixDerivative(N, tau)
    BN_minus1 = bernsteinMatrix(N-1, tau);
    BN_left = [zeros(size(BN_minus1, 1), 1), BN_minus1];
    BN_right = [BN_minus1, zeros(size(BN_minus1, 1), 1)];
    BN_dot = N * (BN_left - BN_right);
end
function FN = FourierMatrix(N, time)
    t = time(:);
    Nt = numel(t);
    T  = time(end);
    k  = (1:N).';                   % N x 1
    ang = 2*pi*(t/T) * k.';         % Nt x N
    FN = ones(Nt, 2*N + 1);
    FN(:, 2:(N+1))     = cos(ang);
    FN(:, (N+2):end)   = sin(ang);
end
function Dm = FourierDiff(N, time)
    t = time(:);
    Nt = numel(t);
    T  = time(end);
    k  = (1:N).';
    ang = 2*pi*(t/T) * k.';         % Nt x N
    Dm = zeros(Nt, 2*N + 1);
    coef = (2*pi/T) * k.';          % 1 x N, broadcasts across Nt
    Dm(:, 2:(N+1))     = -sin(ang) .* coef;
    Dm(:, (N+2):end)   =  cos(ang) .* coef;
end
%% Bernstein differentiation with degree elevation
function Dm = Diff_elev(N, tf)
    Dm_raw = Diff(N, tf);
    Telev = deg_elev(N);
    Dm = Dm_raw * Telev{N-1};
end
function Dm = Diff(N, tf)
    Dm = -[N/tf*eye(N); zeros(1,N)] + [zeros(1,N); N/tf*eye(N)];
end
function Telev = deg_elev(N)
    if N < 5
        error('The approximation order should be at least 5');
    end
    for i = 1:N
        Telev{i} = zeros(i+2, i+1);
        for j = 1:i+1
            Telev{i}(j, j) = i + 1 - (j-1);
            Telev{i}(j+1, j) = (j-1) + 1;
        end
        Telev{i} = (1/(i+1)) * Telev{i}';
    end
end
