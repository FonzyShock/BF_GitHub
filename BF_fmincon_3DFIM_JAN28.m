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

data = struct();
% UAV parameters
data.V_A   = 30;                      % m/s, UAV constant velocity
data.z_A   = 100;                     % m, UAV altitude
data.psi0  = 0;                       % rad, UAV initial heading
data.x0_A  = [0, 0];                  % m, UAV initial position
data.u0    = [cos(data.psi0), sin(data.psi0)];  % initial control (heading direction)

% Target parameters (MOVING TARGET)
data.x_T0  = [400, 0, 0];                % m, target initial position
data.v_T   = [5, 0];                  % m/s, target velocity (known)

% Measurement parameters
data.Sigma = diag([3e-3, 3e-3].^2);   % Measurement noise covariance (rad^2)
data.freq  = 10;                      % Hz, measurement frequency

% Constraint parameters
data.r_NFZ = 100;                     % m, no-fly zone radius
data.kappa = 1e-2;                    % Barrier parameter (tunable)
data.phi_max = deg2rad(30); 
data.g = 9.81;
data.turn_rate_max = (data.g * tan(data.phi_max)) / data.V_A;

% Time parameters
data.tf = 90;                         % s, scenario duration

% Approximation orders
data.N_B = 40;                        % Bernstein polynomial order
data.N_F = 20;                        % Fourier series order

% State and control dimensions
data.Nx = 2;                          % number of states (x, y)
data.Nu = 2;                          % number of controls (cos(psi), sin(psi))

%% ========================================================================
%  PART 1: MIXED BERNSTEIN-FOURIER WITH MOVING TARGET
%  ========================================================================
fprintf('------------------------------------------------------------------\n');
fprintf('  PART 1: Mixed Bernstein-Fourier Optimization (Moving Target)\n');
fprintf('------------------------------------------------------------------\n');

% Total number of coefficients
data.N_tot = data.N_B + 1 + 2*data.N_F + 1;

% Time nodes
data.Nt_BF = data.N_B + 1;
data.time_BF   = linspace(0, 1, data.Nt_BF);

% Construct basis matrices
data.BN_BF  = bernsteinMatrix(data.N_B, data.time_BF);
data.FN_BF  = FourierMatrix(data.N_F, data.time_BF);
data.DBm_BF = bernsteinMatrixDerivative(data.N_B, data.time_BF);
data.DFm_BF = FourierDiff(data.N_F, data.time_BF);

fprintf('Basis matrix condition number: %.2e\n', cond([data.BN_BF, data.FN_BF]));

% Initial guess: fly toward target, then orbit
Cx_BF = zeros(data.N_tot, data.Nx);
Cu_BF = zeros(data.N_tot, data.Nu);

for k = 1:data.Nt_BF
    t_k = data.time_BF(k) * data.tf;
    x_T_k = data.x_T0(1) + data.v_T(1) * t_k;
    y_T_k = data.x_T0(2) + data.v_T(2) * t_k;
    
    t_approach = norm(data.x_T0 + [-data.r_NFZ, 0,0]) / data.V_A;
    
    if t_k <= t_approach
        Cx_BF(k, :) = [data.V_A * t_k, 0];
        Cu_BF(k, :) = [1, 0];
    else
        angle = pi - data.V_A / data.r_NFZ * (t_k - t_approach);
        Cx_BF(k, :) = [x_T_k, y_T_k] + data.r_NFZ * [cos(angle), sin(angle)];
        psi_k = pi/2 - data.V_A / data.r_NFZ * (t_k - t_approach);
        Cu_BF(k, :) = [cos(psi_k), sin(psi_k)];
    end
end

% Reshape for optimizer
X0_BF = [reshape(Cx_BF, [data.Nx * data.N_tot, 1]); reshape(Cu_BF, [data.Nu * data.N_tot, 1])];

% Equality constraints (initial conditions)
% Access BN and FN from data struct
Aeq_BF = blkdiag([data.BN_BF(1,:), data.FN_BF(1,:)], [data.BN_BF(1,:), data.FN_BF(1,:)], ...
                 [data.BN_BF(1,:), data.FN_BF(1,:)], [data.BN_BF(1,:), data.FN_BF(1,:)]);
beq_BF = [data.x0_A'; data.u0'];

% Bounds
lb_BF = -inf(size(X0_BF));
ub_BF =  inf(size(X0_BF));

% Optimizer options
options_BF = optimoptions(@fmincon, ...
    'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1e6, ...
    'OptimalityTolerance', 1e-3,...
    'ConstraintTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'MaxIterations', 5e3, ...
    'Display', 'iter');

% Define cost and constraint functions
costFunc_BF = @(X) costFunc_MixedBF_MovingTarget_3DFIM(X, data);
nonlcon_BF  = @(X) nonlcon_MixedBF_MovingTarget(X, data);

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

% Add Bernstein-specific parameters to the data struct
data.N_bern = 100;
data.r_NFZ_bern = 30; % Relaxed constraint for baseline

% Time nodes
Nt_bern = data.N_bern + 1;
data.time_bern = linspace(0, 1, Nt_bern);

% Basis matrices
data.BN_bern = bernsteinMatrix(data.N_bern, data.time_bern);
data.Dm_bern = Diff_elev(data.N_bern, 1);

% Initial guess
Cx_bern = zeros(data.N_bern + 1, data.Nx);
Cu_bern = zeros(data.N_bern + 1, data.Nu);

for k = 1:Nt_bern
    t_k = data.time_bern(k) * data.tf;
    x_T_k = data.x_T0(1) + data.v_T(1) * t_k;
    y_T_k = data.x_T0(2) + data.v_T(2) * t_k;
    
    t_approach = norm(data.x_T0 + [-data.r_NFZ_bern, 0, 0]) / data.V_A;
    
    if t_k <= t_approach
        Cx_bern(k, :) = [data.V_A * t_k, 0];
        Cu_bern(k, :) = [1, 0];
    else
        angle = pi - data.V_A / data.r_NFZ_bern * (t_k - t_approach);
        Cx_bern(k, :) = [x_T_k, y_T_k] + data.r_NFZ_bern * [cos(angle), sin(angle)];
        psi_k = pi/2 - data.V_A / data.r_NFZ_bern * (t_k - t_approach);
        Cu_bern(k, :) = [cos(psi_k), sin(psi_k)];
    end
end

X0_bern = [reshape(Cx_bern, [data.Nx * (data.N_bern + 1), 1]); ...
           reshape(Cu_bern, [data.Nu * (data.N_bern + 1), 1])];

% Equality constraints
Aeq_bern = zeros(data.Nx + data.Nu, (data.N_bern + 1) * (data.Nx + data.Nu));
Aeq_bern(1, 1) = 1;
Aeq_bern(2, data.N_bern + 2) = 1;
Aeq_bern(3, (data.N_bern + 1) * data.Nx + 1) = 1;
Aeq_bern(4, (data.N_bern + 1) * data.Nx + data.N_bern + 2) = 1;
beq_bern = [data.x0_A'; data.u0'];

% Bounds
lb_bern = -inf(size(X0_bern));
ub_bern =  inf(size(X0_bern));

% Optimizer options
options_bern = optimoptions(@fmincon, ...
    'Algorithm', 'sqp', ...
    'OptimalityTolerance',1e-3, ...
    'MaxFunctionEvaluations', 1e6, ...
    'ConstraintTolerance', 1e-6, ...
    'StepTolerance', 1e-8, ...
    'MaxIterations', 5e3, ...
    'Display', 'iter');

% Define cost and constraint functions
costFunc_bern = @(X) costFunc_Bernstein_MovingTarget_3DFIM(X, data);
nonlcon_bern  = @(X) nonlcon_Bernstein_MovingTarget(X, data);

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

% Evaluation grid based on camera sample frequency
Nt_eval = data.tf * data.freq + 1;
time_eval = linspace(0, 1, Nt_eval);
time_eval_sec = time_eval(:) * data.tf;

% Evaluation basis matrices
% --- For Mixed BF ---
data_eval_BF = data; 
data_eval_BF.BN_BF = bernsteinMatrix(data.N_B, time_eval);
data_eval_BF.FN_BF = FourierMatrix(data.N_F, time_eval);

% --- For Bernstein ---
data_eval_bern = data;
data_eval_bern.BN_bern = bernsteinMatrix(data.N_bern, time_eval);

% 3. Extract Trajectories (Using the modified data structs)
[X_A_BF, Y_A_BF, cos_U_BF, sin_U_BF] = ...
    getFullStates_MixedBF(results.mixedBF.xOut, data_eval_BF);

[X_A_bern, Y_A_bern, cos_U_bern, sin_U_bern] = ...
    getStates_Bernstein(results.bernstein.xOut, data_eval_bern);

% Compute target trajectory

% Compute target trajectory
X_T_traj = data.x_T0(1) + data.v_T(1) * time_eval_sec';
Y_T_traj = data.x_T0(2) + data.v_T(2) * time_eval_sec';

X_A_BF = X_A_BF(:); Y_A_BF = Y_A_BF(:);
X_A_bern = X_A_bern(:); Y_A_bern = Y_A_bern(:);

% Compute distances to target
dist_BF = sqrt((X_A_BF - X_T_traj).^2 + (Y_A_BF - Y_T_traj).^2);
dist_bern = sqrt((X_A_bern - X_T_traj).^2 + (Y_A_bern - Y_T_traj).^2);

% Compute CRLB trace for both methods 
[CRLB_BF, CRLB_BF_cumulative] = computeCRLB_MovingTarget_3DFIM(X_A_BF, Y_A_BF, data.x_T0, data.v_T, data.z_A, data.Sigma, data.tf, time_eval);
[CRLB_bern, CRLB_bern_cumulative] = computeCRLB_MovingTarget_3DFIM(X_A_bern, Y_A_bern, data.x_T0, data.v_T, data.z_A, data.Sigma, data.tf, time_eval);

% --- Calculate FOV Indicators & Bank Angles for Plotting ---
% Mixed BF Bank Angle
psi_BF = unwrap(atan2(sin_U_BF, cos_U_BF));
phi_BF = atan((data.V_A .* gradient(psi_BF, time_eval_sec)) / 9.81);
turn_rate_BF = gradient(psi_BF, time_eval_sec);
fov_BF = get_FOV_Weight(X_A_BF, Y_A_BF, data.z_A, cos_U_BF, sin_U_BF, turn_rate_BF, X_T_traj, Y_T_traj, data.V_A);

% Bernstein Bank Angle
psi_bern = unwrap(atan2(sin_U_bern, cos_U_bern));
phi_bern = atan((data.V_A .* gradient(psi_bern, time_eval_sec)) / 9.81);
turn_rate_bern = gradient(psi_bern, time_eval_sec);
fov_bern = get_FOV_Weight(X_A_bern, Y_A_bern, data.z_A, cos_U_bern, sin_U_bern, turn_rate_bern, X_T_traj, Y_T_traj, data.V_A);

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
results.eval.phi_BF = phi_BF;
results.eval.phi_bern = phi_bern;
results.eval.fov_BF = fov_BF;
results.eval.fov_bern = fov_bern;

fprintf('\n--- Performance Comparison ---\n');
fprintf('%-30s %15s %15s\n', 'Metric', 'Mixed BF', 'Bernstein');
fprintf('%-30s %15s %15s\n', repmat('-',1,30), repmat('-',1,15), repmat('-',1,15));
fprintf('%-30s %15.5e %15.5e\n', 'Final CRLB trace', CRLB_BF, CRLB_bern);
fprintf('%-30s %15.2f %15.2f\n', 'Computation time (s)', results.mixedBF.time_cold, results.bernstein.time);
fprintf('%-30s %15d %15d\n', 'Decision variables', results.mixedBF.numVars, results.bernstein.numVars);

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

% %% ========================================================================
% %  PART 4: COVECTOR MAPPING VALIDATION
% %  ========================================================================
% fprintf('\n------------------------------------------------------------------\n');
% fprintf('  PART 4: Covector Mapping Validation\n');
% fprintf('------------------------------------------------------------------\n');
% 
% Nt_opt = data.Nt_BF;
% if isfield(results.mixedBF.lambda, 'eqnonlin') && ~isempty(results.mixedBF.lambda.eqnonlin)
%     lambda_eq = results.mixedBF.lambda.eqnonlin;
%     lambda_x_NLP = lambda_eq(1:Nt_opt);
%     lambda_y_NLP = lambda_eq(Nt_opt+1:2*Nt_opt);
% 
%     w = 1 / (Nt_opt + 1);
%     lambda_x_scaled = lambda_x_NLP / w;
%     lambda_y_scaled = lambda_y_NLP / w;
% 
%     results.covector.lambda_x_NLP = lambda_x_scaled;
%     results.covector.lambda_y_NLP = lambda_y_scaled;
% else
%     results.covector.lambda_x_NLP = [];
%     results.covector.lambda_y_NLP = [];
% end
% 
% [H_traj, H_components] = computeHamiltonian_MovingTarget_3DFIM(...
%     X_A_BF, Y_A_BF, cos_U_BF, sin_U_BF, ...
%     data_eval_BF.x_T0, data_eval_BF.v_T, data_eval_BF.z_A, data_eval_BF.Sigma, data_eval_BF.tf, data_eval_BF.r_NFZ, data_eval_BF.kappa, data_eval_BF.V_A, ...
%     results.covector.lambda_x_NLP, results.covector.lambda_y_NLP, ...
%     data_eval_BF.time_BF, data_eval_BF.BN_BF, data_eval_BF.FN_BF, data_eval_BF.N_B, data_eval_BF.N_F);
% 
% results.covector.H = H_traj;
% results.covector.H_components = H_components;
% 
% fprintf('\nHamiltonian analysis (low-order, N_B=%d, N_F=%d):\n', data.N_B, data.N_F);
% fprintf('  Mean: %.6f\n', mean(H_traj));
% fprintf('  Std:  %.6f\n', std(H_traj));

%% ========================================================================
%  PART 6: GENERATE FIGURES
%  ========================================================================
fprintf('\n------------------------------------------------------------------\n');
fprintf('  PART 6: Generate Figures\n');
fprintf('------------------------------------------------------------------\n');

% Figure 5: Top-down trajectories
figure('Name', 'TopDown_Trajectories', 'Position', [100, 100, 800, 600]);
hold on; grid on; axis equal;
plot(X_A_BF, Y_A_BF, 'b-', 'LineWidth', 2, 'DisplayName', 'Mixed BF');
plot(X_A_bern, Y_A_bern, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Bernstein-only');
plot(X_T_traj, Y_T_traj, 'r--', 'LineWidth', 2, 'DisplayName', 'Target');
plot(data.x_T0(1), data.x_T0(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Target start');
plot(X_T_traj(end), Y_T_traj(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Target end');

theta_circle = linspace(0, 2*pi, 100);
plot(data.x_T0(1) + data.r_NFZ * cos(theta_circle), data.x_T0(2) + data.r_NFZ * sin(theta_circle), 'r:', 'LineWidth', 1.5);
plot(X_T_traj(end) + data.r_NFZ * cos(theta_circle), Y_T_traj(end) + data.r_NFZ * sin(theta_circle), 'm:', 'LineWidth', 1.5, 'DisplayName', 'NFZ (End)');
plot(data.x0_A(1), data.x0_A(2), 'b^', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'UAV start');
xlabel('x, m'); ylabel('y, m');
legend('Location', 'best');
title(sprintf('Optimal UAV trajectories (3D FIM Trace: BF=%.2e, Bern=%.2e)', CRLB_BF, CRLB_bern));

% Figure 7: Constraint and information
figure('Name', 'Fig7_Constraint_Information', 'Position', [100, 100, 800, 500]);
subplot(2, 1, 1);
hold on; grid on;
plot(time_eval_sec, dist_BF, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Mixed BF');
plot(time_eval_sec, dist_bern, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Bernstein-only');
yline(data.r_NFZ, 'r--', 'NFZ');
xlabel('Time, s'); ylabel('Distance to target, m');
legend('Location', 'best');
title('UAV-Target distance');

subplot(2, 1, 2);
hold on; grid on;
plot(time_eval_sec, CRLB_BF_cumulative, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Mixed BF');
plot(time_eval_sec, CRLB_bern_cumulative, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Bernstein-only');
xlabel('Time, s'); ylabel('Trace(FIM^{-1})');
legend('Location', 'best');
title('Cumulative A-Optimality Cost');

% New Figure: State and Control History
figure('Name', 'State_Control_History', 'Position', [150, 150, 800, 800]);

% X Position
subplot(4, 1, 1);
hold on; grid on;
plot(time_eval_sec, X_A_BF, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Mixed BF');
plot(time_eval_sec, X_A_bern, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Bernstein');
ylabel('X [m]');
title('UAV State and Control History');
legend('Location', 'best');

% Y Position
subplot(4, 1, 2);
hold on; grid on;
plot(time_eval_sec, Y_A_BF, 'b-', 'LineWidth', 1.5);
plot(time_eval_sec, Y_A_bern, 'g--', 'LineWidth', 1.5);
ylabel('Y [m]');

% Heading (Psi)
subplot(4, 1, 3);
hold on; grid on;
plot(time_eval_sec, rad2deg(psi_BF), 'b-', 'LineWidth', 1.5);
plot(time_eval_sec, rad2deg(psi_bern), 'g--', 'LineWidth', 1.5);
ylabel('Heading \psi [deg]');

% Bank Angle (Phi)
subplot(4, 1, 4);
hold on; grid on;
plot(time_eval_sec, rad2deg(phi_BF), 'b-', 'LineWidth', 1.5);
plot(time_eval_sec, rad2deg(phi_bern), 'g--', 'LineWidth', 1.5);
yline(30, 'r:', 'Limit'); yline(-30, 'r:');
ylabel('Bank Angle \phi [deg]');
xlabel('Time [s]');

% FOV Visibility Check
figure('Name', 'FOV_Check', 'Position', [200, 200, 800, 400]);
hold on; grid on;
plot(time_eval_sec, fov_BF, 'b-', 'LineWidth', 2, 'DisplayName', 'Mixed BF');
plot(time_eval_sec, fov_bern, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Bernstein');
title('Target Visibility (FOV Weight)');
xlabel('Time (s)'); ylabel('Weight (0-1)');
legend; ylim([-0.1 1.1]);

% Save results
save('moving_target_results.mat', 'results');

%% ========================================================================
%  DYNAMICS FEASIBILITY CHECK
%  ========================================================================
checkFeasibility('Mixed BF', X_A_BF, Y_A_BF, cos_U_BF, sin_U_BF, ...
                 time_eval, data.V_A, data.tf);

checkFeasibility('Bernstein', X_A_bern, Y_A_bern, cos_U_bern, sin_U_bern, ...
                 time_eval, data.V_A, data.tf);
%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

%% Cost function for Mixed BF (3D FIM Accumulated)
function J = costFunc_MixedBF_MovingTarget_3DFIM(X, data)
    N_B = data.N_B; N_F = data.N_F; N_tot = data.N_tot;
    Nx = data.Nx; Nu = data.Nu;
    BN = data.BN_BF; FN = data.FN_BF; 
    DBm = data.DBm_BF; DFm = data.DFm_BF;
    x_T0 = data.x_T0; v_T = data.v_T;
    z_A = data.z_A; V_A = data.V_A;
    Sigma = data.Sigma; tf = data.tf;
    r_NFZ = data.r_NFZ; kappa = data.kappa;
    time = data.time_BF;
   
    % Extract Coefficients
    Cx = reshape(X(1:Nx*N_tot), [N_tot, Nx]);
    Cu = reshape(X(Nx*N_tot+1:end), [N_tot, Nu]);
    
    CxB = Cx(1:N_B+1, :); CxF = Cx(N_B+2:end, :);
    CuB = Cu(1:N_B+1, :); CuF = Cu(N_B+2:end, :);
    
    % Reconstruct states and controls 
    X_A = BN * CxB(:,1) + FN * CxF(:,1);
    Y_A = BN * CxB(:,2) + FN * CxF(:,2);
    cos_U = BN * CuB(:,1) + FN * CuF(:,1);
    sin_U = BN * CuB(:,2) + FN * CuF(:,2);

    % Calculate Turn Rate for FOV Weighting
    d_cos = DBm * CuB(:,1) + DFm * CuF(:,1);
    d_sin = DBm * CuB(:,2) + DFm * CuF(:,2);
    turn_rate = (cos_U .* d_sin - sin_U .* d_cos) / tf;
    
    % Target Trajectory
    tsec = (time(:) * tf);
    X_T  = x_T0(1) + v_T(1) * tsec;
    Y_T  = x_T0(2) + v_T(2) * tsec;
    Z_T = x_T0(3); 
    
    % Calculate FOV Weights
    %FOV_Ind = get_FOV_Weight(X_A, Y_A, z_A, cos_U, sin_U, turn_rate, X_T, Y_T, V_A);
FOV_Ind = 1;
    % 3D Relative Geometry (r = Target - Cam)
    rx = X_T - X_A;
    ry = Y_T - Y_A;
    rz = Z_T - z_A; % This will be negative altitude
    
    % Relative Geometry
    tol = 1e-6;
    rxy2 = rx.^2 + ry.^2 + tol;
    rxy  = sqrt(rxy2);
    r2   = rxy2 + rz.^2 + tol;
    
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
    
    % Accumulate FIM (Weighted)
    F11 = sum(FOV_Ind .* (w_az * h_az_x.^2 + w_el * h_el_x.^2));
    F22 = sum(FOV_Ind .* (w_az * h_az_y.^2 + w_el * h_el_y.^2));
    F33 = sum(FOV_Ind .* (w_az * h_az_z.^2 + w_el * h_el_z.^2));
    F12 = sum(FOV_Ind .* (w_az * h_az_x .* h_az_y + w_el * h_el_x .* h_el_y));
    F13 = sum(FOV_Ind .* (w_az * h_az_x .* h_az_z + w_el * h_el_x .* h_el_z));
    F23 = sum(FOV_Ind .* (w_az * h_az_y .* h_az_z + w_el * h_el_y .* h_el_z));
    
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
function J = costFunc_Bernstein_MovingTarget_3DFIM(X, data)
    N = data.N_bern; 
    Nx = data.Nx; 
    Nu = data.Nu;
    BN = data.BN_bern; 
    Dm = data.Dm_bern;
    x_T0 = data.x_T0; 
    v_T = data.v_T;
    z_A = data.z_A; 
    V_A = data.V_A;
    Sigma = data.Sigma; 
    tf = data.tf;
    r_NFZ = data.r_NFZ_bern;
    kappa = data.kappa;
    time = data.time_bern;

    % Extract Coefficients
    Cx = reshape(X(1:Nx*(N+1)), [N+1, Nx]);
    Cu = reshape(X(Nx*(N+1)+1:end), [N+1, Nu]);
    
    % Reconstruct
    X_A = BN * Cx(:, 1);
    Y_A = BN * Cx(:, 2);
    cos_U = BN * Cu(:, 1);
    sin_U = BN * Cu(:, 2);
    
    % Turn Rate
    d_cos = Dm' * Cu(:, 1);
    d_sin = Dm' * Cu(:, 2);
    turn_rate = (cos_U .* d_sin - sin_U .* d_cos) / tf;

    % Target Trajectory
    tsec = time(:) * tf;
    X_T = x_T0(1) + v_T(1) * tsec;
    Y_T = x_T0(2) + v_T(2) * tsec;
    Z_T = x_T0(3);
    
    % Calculate FOV Weights
    %FOV_Ind = get_FOV_Weight(X_A, Y_A, z_A, cos_U, sin_U, turn_rate, X_T, Y_T, V_A);
    FOV_Ind = 1;
    % 3D Relative Geometry (r = Target - Cam)
    rx = X_T - X_A;
    ry = Y_T - Y_A;
    rz = Z_T - z_A; % This will be negative altitude
    
    % Relative Geometry
    tol = 1e-6;
    rxy2 = rx.^2 + ry.^2 + tol;
    rxy  = sqrt(rxy2);
    r2   = rxy2 + rz.^2 + tol;
    
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
    
    % Accumulate FIM (Weighted)
    F11 = sum(FOV_Ind .* (w_az * h_az_x.^2 + w_el * h_el_x.^2));
    F22 = sum(FOV_Ind .* (w_az * h_az_y.^2 + w_el * h_el_y.^2));
    F33 = sum(FOV_Ind .* (w_az * h_az_z.^2 + w_el * h_el_z.^2));
    F12 = sum(FOV_Ind .* (w_az * h_az_x .* h_az_y + w_el * h_el_x .* h_el_y));
    F13 = sum(FOV_Ind .* (w_az * h_az_x .* h_az_z + w_el * h_el_x .* h_el_z));
    F23 = sum(FOV_Ind .* (w_az * h_az_y .* h_az_z + w_el * h_el_y .* h_el_z));
    
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

%% Compute CRLB History (3D FIM Accumulated)
function [CRLB_final, CRLB_cumulative] = computeCRLB_MovingTarget_3DFIM(X_A, Y_A, x_T0, v_T, z_A, Sigma, tf, time)
    
    Nt = length(time);
    tsec = time * tf;
    
    X_T = x_T0(1) + v_T(1) * tsec(:);
    Y_T = x_T0(2) + v_T(2) * tsec(:);
    Z_T = x_T0(3);
    
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
       
    CRLB_cumulative = zeros(Nt, 1);
    FIM = zeros(3);
    
    for k = 1:Nt
        J_k = [h_az_x(k), h_az_y(k), h_az_z(k); 
               h_el_x(k), h_el_y(k), h_el_z(k)];
        FIM = FIM + J_k' * Sigma * J_k;
            
        if k < 5 || rcond(FIM) < 1e-12
            CRLB_cumulative(k) = NaN; % Not observable yet
        else
            CRLB_cumulative(k) = trace(inv(FIM));
        end
    end
    CRLB_final = CRLB_cumulative(end);
end

%% FOV Weight Calculator (Analytic Bank Angle)
function weights = get_FOV_Weight(X_A, Y_A, z_A, cos_U, sin_U, turn_rate, X_T, Y_T, V_A)
    % Parameters
    g = 9.81;
    lim_az    = deg2rad(60); 
    lim_el_up = deg2rad(45);
    lim_el_dn = deg2rad(-45);
    k = 10; EXP_MIN = -700; EXP_MAX = 700;

    % Bank Angle
    phi = atan( (V_A .* turn_rate) / g );
    cPhi = cos(phi); sPhi = sin(phi);

    % Geometry
    dx = X_T - X_A;
    dy = Y_T - Y_A;
    dz = z_A; 

    % Rotation (Yaw)
    x1 =  cos_U .* dx + sin_U .* dy;
    y1 = -sin_U .* dx + cos_U .* dy;
    z1 =  dz;

    % Rotation (Roll)
    Xb = x1;
    Yb =  cPhi .* y1 + sPhi .* z1;
    Zb = -sPhi .* y1 + cPhi .* z1;

    % Angles
    rel_az = atan2(Yb, Xb);
    range_xy = sqrt(Xb.^2 + Yb.^2);
    rel_el = atan2(-Zb, range_xy);

    % Sigmoids
    clip = @(v) min(max(v, EXP_MIN), EXP_MAX);
    sig  = @(x, thr) 1 ./ (1 + exp(clip(k * (x - thr))));

    w_az = sig(abs(rel_az), lim_az/2);
    w_el = sig(rel_el, lim_el_up) .* (1 - sig(rel_el, lim_el_dn));

    weights = w_az .* w_el;
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
function [c, ceq] = nonlcon_MixedBF_MovingTarget(X, data)

    N_B = data.N_B;
    N_F = data.N_F;
    N_tot = data.N_tot;
    Nx  = data.Nx; 
    Nu  = data.Nu;
    V_A = data.V_A;
    tf  = data.tf;
    DBm = data.DBm_BF; 
    DFm = data.DFm_BF;
    BN = data.BN_BF; 
    FN = data.FN_BF;
    turn_rate_max = data.turn_rate_max;
    
    % Extract coefficients
    Cx = reshape(X(1:Nx*N_tot), [N_tot, Nx]);
    Cu = reshape(X(Nx*N_tot+1:end), [N_tot, Nu]);
    
    % Controls
    cos_U = BN * Cu(1:N_B+1, 1) + FN * Cu(N_B+2:end, 1);
    sin_U = BN * Cu(1:N_B+1, 2) + FN * Cu(N_B+2:end, 2);
    
    % State derivatives
    dX_A = DBm * Cx(1:N_B+1, 1) + DFm * Cx(N_B+2:end, 1);
    dY_A = DBm * Cx(1:N_B+1, 2) + DFm * Cx(N_B+2:end, 2);

    % Turn Rate
    d_cos = DBm * Cu(1:N_B+1, 1) + DFm * Cu(N_B+2:end, 1);
    d_sin = DBm * Cu(1:N_B+1, 2) + DFm * Cu(N_B+2:end, 2);
    turn_rate = (cos_U .* d_sin - sin_U .* d_cos) / tf;
    
    % Dynamics constraints
    ceq = [dX_A - tf * V_A * cos_U; 
           dY_A - tf * V_A * sin_U; 
           cos_U.^2 + sin_U.^2 - 1];
    
    c = abs(turn_rate) - turn_rate_max;
end

%% Nonlinear constraints for Bernstein
function [c, ceq] = nonlcon_Bernstein_MovingTarget(X, data)

    N = data.N_bern; Nx = data.Nx; Nu = data.Nu;
    V_A = data.V_A; tf = data.tf;
    Dm = data.Dm_bern;
    turn_rate_max = data.turn_rate_max;

    % Extract coefficients
    Cx = reshape(X(1:Nx*(N+1)), [N+1, Nx]);
    Cu = reshape(X(Nx*(N+1)+1:end), [N+1, Nu]);
    
    % Controls
    cos_U = Cu(:, 1);
    sin_U = Cu(:, 2);
    
    % State derivatives
    dX_A = Dm' * Cx(:, 1);
    dY_A = Dm' * Cx(:, 2);

    % Turn Rate
    d_cos = Dm' * cos_U;
    d_sin = Dm' * sin_U;
    turn_rate = (cos_U .* d_sin - sin_U .* d_cos) / tf;
    
    ceq = [dX_A - tf * V_A * cos_U; 
           dY_A - tf * V_A * sin_U; 
           cos_U.^2 + sin_U.^2 - 1];

    c = abs(turn_rate) - turn_rate_max;
end

%% Dynamics Feasibility Check
function checkFeasibility(methodName, X, Y, cos_U, sin_U, time, V_A, tf)
    dx = gradient(X);
    dy = gradient(Y);
    dist_traj = sum(sqrt(dx.^2 + dy.^2));
    dist_expected = V_A * tf;
    
    x_integrated = X(1) + cumtrapz(time * tf, V_A * cos_U);
    y_integrated = Y(1) + cumtrapz(time * tf, V_A * sin_U);
    
    fprintf('\n--- Feasibility Check: %s ---\n', methodName);
    fprintf('  Expected Dist: %.2f m | Actual: %.2f m (Err: %.4f%%)\n', ...
        dist_expected, dist_traj, abs(dist_traj - dist_expected)/dist_expected*100);
    fprintf('  Integ. Err (X,Y): (%.2e, %.2e)\n', ...
        abs(X(end) - x_integrated(end)), abs(Y(end) - y_integrated(end)));
end

%% Get full states for Mixed BF
function [X_A, Y_A, cos_U, sin_U] = getFullStates_MixedBF(X, data)
    N_B = data.N_B;
    N_tot = data.N_tot;
    Nx = data.Nx;
    BN = data.BN_BF;
    FN = data.FN_BF;
    Nu = data.Nu;

    Cx = reshape(X(1:Nx*N_tot), [N_tot, Nx]);
    Cu = reshape(X(Nx*N_tot+1:end), [N_tot, Nu]);
    
    X_A = BN * Cx(1:N_B+1, 1) + FN * Cx(N_B+2:end, 1);
    Y_A = BN * Cx(1:N_B+1, 2) + FN * Cx(N_B+2:end, 2);
    cos_U = BN * Cu(1:N_B+1, 1) + FN * Cu(N_B+2:end, 1);
    sin_U = BN * Cu(1:N_B+1, 2) + FN * Cu(N_B+2:end, 2);
end

%% Get states and controls for Bernstein
function [X_A, Y_A, cos_U, sin_U] = getStates_Bernstein(X, data)
    N = data.N_bern;
    Nx = data.Nx;
    Nu = data.Nu;
    BN = data.BN_bern;

    Cx = reshape(X(1:Nx*(N+1)), [N+1, Nx]);
    Cu = reshape(X(Nx*(N+1)+1:end), [N+1, Nu]);
    X_A = BN * Cx(:, 1); 
    Y_A = BN * Cx(:, 2);
    cos_U = BN * Cu(:, 1); 
    sin_U = BN * Cu(:, 2);
end

%% Basis Functions

% Bernstein matrix
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

% Bernstein derivative matrix
function BN_dot = bernsteinMatrixDerivative(N, tau)
    BN_minus1 = bernsteinMatrix(N-1, tau);
    BN_left = [zeros(size(BN_minus1, 1), 1), BN_minus1];
    BN_right = [BN_minus1, zeros(size(BN_minus1, 1), 1)];
    BN_dot = N * (BN_left - BN_right);
end

% Fourier matrix
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

% Fourier differentiation matrix
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
