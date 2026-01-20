% Federated_EKF
% Bearing-only federated EKF for supporting UAVs with 3D stationary target state [x y z].
% - Uses results_Coop/data_Coop from optimization
% - Outputs RMSE struct with boat/target fields into workspace

%%%%%% Load Data %%%%%%%%%%% 
if ~(exist('data_Coop','var') && exist('results_Coop','var'))
    [~, computerName] = system('hostname');
    computerName = strtrim(computerName);
    loadFiles = struct();
    switch computerName
        case 'DESKTOP'
            loadFiles.folderPath = 'C:\GitRepos\COMET\PhD Research\Output Files\250711';
        case 'Surface4'
            loadFiles.folderPath = 'C:\Users\alfon\OneDrive\Documents\MATLAB\Paper1\Output Files\250311';
        otherwise
            warning('Unknown host and no in-memory data; add path or pre-load data_Coop/results_Coop.');
            loadFiles = [];
    end
    if ~isempty(loadFiles)
        fileNames = {'data_Coop_5UAV.mat','results_Coop_5UAV.mat'};
        for kf = 1:numel(fileNames)
            f = fullfile(loadFiles.folderPath, fileNames{kf});
            if exist(f,'file')
                S = load(f);
                vars = fieldnames(S);
                for v = 1:numel(vars)
                    assignin('caller', vars{v}, S.(vars{v}));
                end
                fprintf('Loaded %s\n', fileNames{kf});
            else
                warning('File not found: %s', f);
            end
        end
    end
else
    fprintf('Using in-memory data_Coop/results_Coop.\n');
end

%%%%%%% Setup %%%%%%%%%%
data        = data_Coop;
results     = results_Coop;
num_uav     = data.num_uav;
num_targets = data.num_targets;
dt          = 1 / data.sample_rate;
alpha_ipopt = results.xout(end);
tf          = (1+alpha_ipopt) * data.tf;

% Labels
if num_targets == 2
    names = {'Support USV','Target USV'};
else
    names = arrayfun(@(i)sprintf('Target %d',i), 1:num_targets, 'uni', 0);
end

%%%%%%%%% Number of Monte Carlo Sims %%%%%%%%%%%%
num_simulations = 1000;

%%%%%%%% EKF model %%%%%%%%%%%%%
I3      = eye(3);
F       = I3;
Q       = data.Q;


sigma_xy0 = 100;
sigma_z0  = 1;  % z is different since there is no change expected. Can be swell height.
P0_fallback        = diag([sigma_xy0^2, sigma_xy0^2, sigma_z0^2]);

% Get prior from main code otherwise fallback
if isfield(data,'P0_AT') && ~isempty(data.P0_AT)
    P0_AT = data.P0_AT;
else
    P0_AT = P0_fallback;
end
if isfield(data,'P0_AB') && ~isempty(data.P0_AB)
    P0_AB = data.P0_AB;
else
    P0_AB = P0_fallback;
end

% Single lambda for all targets
if isfield(data,'lambda') && ~isempty(data.lambda)
    lambda = data.lambda;
else
    lambda = 1;
end
lambda_eps = 1e-12;
epsilon = 1e-9;

%%%%%%%%% Extract cooperative trajectories %%%%%%%%%%%%
results = extractTrajectories_Cooperative(results.xout, data);

UAV_positions = struct();
UAV_yaw       = struct();
UAV_roll      = struct();
traj          = struct();

for idxUAV = 1:num_uav
    Xi  = results.(sprintf('UAV%d', idxUAV)).X_A;
    Yi  = results.(sprintf('UAV%d', idxUAV)).Y_A;
    Psi = results.(sprintf('UAV%d', idxUAV)).Psi_A;
    Phi = atan(results.(sprintf('UAV%d', idxUAV)).tan_Phi_A);

    UAV_positions(idxUAV).trajectory = [Xi, Yi];
    UAV_yaw(idxUAV).yaw_angle        = Psi;
    UAV_roll(idxUAV).roll_angle      = Phi;

    if isfield(data,'isGimballed') && data.isGimballed
         cPsi_g_A = results.(sprintf('UAV%d', idxUAV)).cPsi_g_A;
         sPsi_g_A = results.(sprintf('UAV%d', idxUAV)).sPsi_g_A;
         cTheta_g_A = results.(sprintf('UAV%d', idxUAV)).cTheta_g_A;
         sTheta_g_A = results.(sprintf('UAV%d', idxUAV)).sTheta_g_A;
         u_Psi_g_A = results.(sprintf('UAV%d', idxUAV)).u_Psi_g_A;
         u_Theta_g_A = results.(sprintf('UAV%d', idxUAV)).u_Theta_g_A;
         traj.(sprintf('UAV%d', idxUAV)) = struct( ...
        'Psi_A', Psi, ... 
        'Phi_A',  Phi, ...
        'X_A',    Xi, ...
        'Y_A',    Yi, ...
        'Z_A',    data.(sprintf('z_A%d', idxUAV)) * ones(size(Xi)), ...
        'cPsi_g_A', cPsi_g_A, ...
        'sPsi_g_A', sPsi_g_A, ...
        'cTheta_g_A', cTheta_g_A, ...
        'sTheta_g_A', sTheta_g_A, ...
        'u_Psi_g_A', u_Psi_g_A, ...
        'u_Theta_g_A', u_Theta_g_A);

    else
    traj.(sprintf('UAV%d', idxUAV)) = struct( ...
        'Psi_A', Psi, ... 
        'Phi_A',  Phi, ...
        'X_A',    Xi, ...
        'Y_A',    Yi, ...
        'Z_A',    data.(sprintf('z_A%d', idxUAV)) * ones(size(Xi)));

    
    end
end

%Min trajectory length across UAVs
T_traj = min(arrayfun(@(i) size(UAV_positions(i).trajectory,1), 1:num_uav));
data.num_time_steps = T_traj;

% %% Normalize to shortest UAV time
% true_target_pos_full = normalize_target_intervals(data);
% assert(size(true_target_pos_full,2) >= num_targets, ...
%     'Truth has %d targets but num_targets=%d.', size(true_target_pos_full,2), num_targets);
%
% true_target_pos = true_target_pos_full(:, 1:num_targets, :);   % [2 x num_targets x Ttruth]
% T_truth = size(true_target_pos, 3);
%
% T   = min(T_traj, T_truth);
% dtv = (0:T-1) * (tf/T);
%
% % num_time_steps      = T;
% % data.num_time_steps = T;
%
% %% Set true target positions and init camera
% true_traj = zeros(3, T, num_targets); % [x;y;z] with z=0
% for tgt_idx = 1:num_targets
%     true_traj(1,:,tgt_idx) = squeeze(true_target_pos(1,tgt_idx,1:T));
%     true_traj(2,:,tgt_idx) = squeeze(true_target_pos(2,tgt_idx,1:T));
%     % z=0
% end

%%%%%%%%%%%% Normalize to shortest UAV time %%%%%%%%%%%
raw = normalize_target_intervals(data);     % Ntgt × 3 × T_raw

% Convert to [2 × Ntgt × T] (keep x,y, drop z)
true_target_pos = permute(raw(:,1:2,:), [2 1 3]);   % 2 × Ntgt × T_raw
Ntgt            = size(true_target_pos, 2);
Ttruth          = size(true_target_pos, 3);

assert(Ntgt >= num_targets, 'Truth has %d targets but num_targets=%d.', Ntgt, num_targets);

T               = min(T_traj, Ttruth);
dtv             = (0:T-1) * (tf/T);
num_time_steps  = T;

%%%%% Build true target trajectories %%%%%%%%%%
true_traj           = zeros(3, T, num_targets); % [3 × T × Ntgt]
true_traj(1:2,:,:)  = permute(true_target_pos(:,1:num_targets,1:T), [1 3 2]);
init_xy             = squeeze(true_target_pos(:,1:num_targets,1)).';  % [Ntgt × 2]
cam_pos_all = cell(1,num_uav); % 3×T each

for uav_idx = 1:num_uav
    X_A = UAV_positions(uav_idx).trajectory(1:T,1).';
    Y_A = UAV_positions(uav_idx).trajectory(1:T,2).';
    Z_A = data.(sprintf('z_A%d',uav_idx)) * ones(1,T);
    cam_pos_all{uav_idx} = [X_A; Y_A; Z_A];
end

%%%%%%%%%%% Monte Carlo Setup %%%%%%%%%%
monte_carlo_results        = struct([]);
monte_carlo_RMSE_over_time = struct([]);
residuals_sq               = zeros(num_targets, T);
%init_xy = reshape(true_target_pos(1:2, :, 1), 2, []).';

if size(init_xy,1) ~= num_targets
    error('Size mismatch: init_xy has %d rows but num_targets=%d.', size(init_xy,1), num_targets);
end

% Per-target cueing errors for x and y
cue_boat   = data.cueing.boat;    % 1-sigma [m] for boat x and y
cue_target = data.cueing.target;  % 1-sigma [m] for target x and y

%%  %%%%%%%% Monte Carlo loop %%%%%%%%%%%%
for sim_idx = 1:num_simulations

    % Noisy initial estimates (independent noise on x and y for each track)
    % init_xy is 2×2: row 1 boat [x y], row 2 target [x y]
    cue_mat = [cue_boat,  cue_boat;    % row 1 applies to boat x,y
               cue_target, cue_target];% row 2 applies to target x,y
    noise_xy = cue_mat .* randn(size(init_xy));
    initial_estimates = [init_xy + noise_xy, zeros(2,1)]; % z initialized to 0

    % Target-specific prior covariance: P0_by_tgt{j} = P0_j / lambda
    P0_boat   = P0_AB / max(lambda, lambda_eps); % track 1
    P0_target = P0_AT / max(lambda, lambda_eps);
    % P0_by_tgt = cell(num_targets,1);


    % Initialize per-UAV EKFs & federated (fed for later code builds)
    EKF_results = repmat(struct('estimates',[],'covariances',[]), num_uav, num_targets);
    for uav_idx = 1:num_uav
        for tgt_idx = 1:num_targets
            EKF_results(uav_idx, tgt_idx).estimates          = zeros(T,3);
            EKF_results(uav_idx, tgt_idx).covariances        = zeros(3,3,T);
            EKF_results(uav_idx, tgt_idx).estimates(1,:)     = initial_estimates(tgt_idx,:);
            EKF_results(uav_idx, tgt_idx).covariances(:,:,1) = iff(tgt_idx==1, P0_boat, P0_target);
        end
    end

    federated_results = repmat(struct('estimates',[],'covariances',[]), 1, num_targets);
    for tgt_idx = 1:num_targets
        federated_results(tgt_idx).estimates                = zeros(T,3);
        federated_results(tgt_idx).covariances              = zeros(3,3,T);
        federated_results(tgt_idx).estimates(1,:)           = initial_estimates(tgt_idx,:);
        federated_results(tgt_idx).covariances(:,:,1)       = iff(tgt_idx==1, P0_boat, P0_target);
    end


    %%%%%%% time loop %%%%%%%%%
    for t = 1:num_time_steps
        for tgt_idx = 1:num_targets
            valid_meas          = false(1, num_uav); % Reset. Used by EKF

            for uav_idx = 1:num_uav
                % Prior
                if t == 1
                    x_est       = EKF_results(uav_idx, tgt_idx).estimates(t,:).';
                    P           = EKF_results(uav_idx, tgt_idx).covariances(:,:,t);
                else
                    x_est       = EKF_results(uav_idx, tgt_idx).estimates(t-1,:).';
                    P           = EKF_results(uav_idx, tgt_idx).covariances(:,:,t-1);
                end

                % Predict
                x_pred          = F*x_est;
                P_pred          = F*P*F' + Q;

                % UAV position at time t
                x_A             = cam_pos_all{uav_idx}(1,t);
                y_A             = cam_pos_all{uav_idx}(2,t);
                z_A             = cam_pos_all{uav_idx}(3,t);

                % True target position (for measurements)
                x_T             = true_traj(1,t,tgt_idx);
                y_T             = true_traj(2,t,tgt_idx);
                z_T             = 0;
                bearing_true    = wrapToPi(atan2(y_T-y_A,x_T - x_A));
                loc_T           = [x_T; y_T; z_T];

                % Expected bearing from the predicted state
                expected_bearing = wrapToPi( atan2( x_pred(2)-y_A, x_pred(1)-x_A));

                % FOV check (All UAVS and t)
                InFOV_val        = InFOV(cam_pos_all{uav_idx}, loc_T, traj, num2str(uav_idx), data);

                % Threshold FOV
                in_fov_t         = InFOV_val(t) > 0.5;%1e-6;

                if in_fov_t
                    % We have a valid measurement for the time step
                    valid_meas(uav_idx) = true;

                    % Linearization
                    dx          = x_pred(1) - x_A;
                    dy          = x_pred(2) - y_A;
                    dz          = x_pred(3) - z_A;

                    rxy2        = dx^2 + dy^2;          
                    rxy         = max(sqrt(rxy2),epsilon);
                    r2          = rxy2 + dz^2;           

                    % Expected measurement
                    az_hat      = atan2(dy, dx);
                    el_hat      = atan2(dz, rxy);
                    z_hat       = [az_hat; el_hat];
                    
                    % True measurement
                    dx_true     = x_T - x_A;
                    dy_true     = y_T - y_A;
                    dz_true     = z_T - z_A;  % z_T = 0 
                    rxy_true    = sqrt(dx_true^2 + dy_true^2);
                    az_true     = wrapToPi(atan2(dy_true, dx_true));
                    el_true     = atan2(dz_true, rxy_true);

                    % Noise (same R as used in update below)
                    R_i         = (tgt_idx-1)*2 + (1:2);
                    Rm          = data.Sigma(R_i, R_i);
                    L           = chol((Rm+Rm')/2, 'lower');
                    v           = L*randn(2,1);

                    z_meas      = [wrapToPi(az_true + v(1));
                                   el_true          + v(2)];

                    % Innovation (wrapped azimuth)
                    residual    = [wrapToPi(z_meas(1) - z_hat(1));
                                   z_meas(2)          - z_hat(2)];

                    % Jacobian H wrt [x y z]
                    H           = zeros(2,3);

                    % Azimuth (theta)
                    H(1,1)      = -dy / max(rxy2, epsilon); % dtheta/dx
                    H(1,2)      =  dx / max(rxy2, epsilon); % dtheta/dy
                    H(1,3)      =  0; % dtheta/dz

                    % Elevation (phi)
                    H(2,1)      = -(dx*dz) / max(rxy*r2, epsilon); % dphi/dx
                    H(2,2)      = -(dy*dz) / max(rxy*r2, epsilon); % dphi/dy
                    H(2,3)      =  rxy     / max(r2, epsilon); % dphi/dz

                    % EKF update (Joseph form for numerical stability)
                    S           = H*P_pred*H' + Rm;
                    K           = P_pred*H'/S;
                    x_upd       = x_pred + K*residual;
                    P_upd       = (I3-K*H)*P_pred*(I3-K*H)' + K*Rm*K';

                    % Store update
                    EKF_results(uav_idx,tgt_idx).estimates(t,:)     = x_upd.';
                    EKF_results(uav_idx,tgt_idx).covariances(:,:,t) = P_upd;
                else
                    % propagate if no measurement
                    EKF_results(uav_idx,tgt_idx).estimates(t,:)     = x_pred.';
                    EKF_results(uav_idx,tgt_idx).covariances(:,:,t) = P_pred;
                end
            end

            %  Terminal-only: NO fusion inside time loop
            if t > 1
                federated_results(tgt_idx).estimates(t,:)     = federated_results(tgt_idx).estimates(t-1,:);
                federated_results(tgt_idx).covariances(:,:,t) = federated_results(tgt_idx).covariances(:,:,t-1);
            else
                federated_results(tgt_idx).estimates(t,:)     = initial_estimates(tgt_idx,:);
                federated_results(tgt_idx).covariances(:,:,t) = iff(tgt_idx==1, P0_boat, P0_target);
            end
        end
    end % End of single run

%%%%%% TERMINAL FUSION with common-info subtraction %%%%%
% Builds a fused posterior at t = T that avoids double-counting the shared prior/dynamics.
% Works per target; keeps earlier times unchanged.

for tgt_idx = 1:num_targets
    
    P0_base  = iff(tgt_idx==1, P0_AB, P0_AT);
    %I_common  = max(lambda, lambda_eps) * inv((P0_base + P0_base')/2); %#ok<MINV>

    P0_local  = P0_base / max(lambda, lambda_eps);
    I_common  = inv((P0_local + P0_local')/2);


    % Prior mean at T for F=I (no control) = sim's initial guess
    x_common            = initial_estimates(tgt_idx,:).';   % 3x1

    % Sum local posteriors in information form
    I_sum               = zeros(3);
    y_sum               = zeros(3,1); % "information-state" = I * x
    
    for uav_idx = 1:num_uav
        PiT             = EKF_results(uav_idx, tgt_idx).covariances(:,:,T);
        xiT             = EKF_results(uav_idx, tgt_idx).estimates(T,:).';
        PiT             = (PiT + PiT.')/2; % symmetrize for inv safety
        Ii              = inv(PiT);
        I_sum           = I_sum + Ii;
        y_sum           = y_sum + Ii*xiT;
    end

    % Subtract duplicated common info that appears in every local posterior:
    % I_fused           = Σ I_post(i) - (M-1) * I_common
    I_fused             = I_sum - (num_uav-1)*I_common;
    I_fused             = (I_fused + I_fused.')/2;  % keep symmetric

    % SPD guard
    [R,p]               = chol(I_fused);
    if p > 0
        % tiny jitter proportional to scale
        s               = max(1, norm(I_fused,'fro')/3);
        I_fused         = I_fused + epsilon*s*I3;
        I_fused         = (I_fused + I_fused.')/2;
    end

    P_fused_T           = inv(I_fused);

    % Fused mean
    % x_fused           = P_fused * ( Σ I_post(i) x_i  − (M−1) I_common x_common )
    y_fused             = y_sum - (num_uav-1)*I_common*x_common;
    x_fused_T           = P_fused_T * y_fused; %#ok<MINV>

    % Terminal fusion result
    federated_results(tgt_idx).estimates(T,:)       = x_fused_T.';
    federated_results(tgt_idx).covariances(:,:,T)   = P_fused_T;
end

   %%% Build per-time RMSE for current sim 
    RMSE_FED_over_time_tmp                          = zeros(num_targets, T);
    for tgt_idx = 1:num_targets
        for t = 1:T
            % gt_xy                                   = [true_traj(1,t,tgt_idx), true_traj(2,t,tgt_idx)];
            % fed_xy                                  = federated_results(tgt_idx).estimates(t,1:2);
            % RMSE_FED_over_time_tmp(tgt_idx,t)       = norm(gt_xy - fed_xy);
            gt_xyz  = [true_traj(1,t,tgt_idx), true_traj(2,t,tgt_idx), 0];
            fed_xyz = federated_results(tgt_idx).estimates(t,1:3);
            RMSE_FED_over_time_tmp(tgt_idx,t) = norm(gt_xyz - fed_xyz);
        end
    end

    % Save current simulation’s results
    monte_carlo_results(sim_idx).EKF_results        = EKF_results;
    monte_carlo_results(sim_idx).federated_results  = federated_results;
    monte_carlo_results(sim_idx).RMSE_overall       = sqrt(mean(RMSE_FED_over_time_tmp.^2,2));
    monte_carlo_RMSE_over_time(sim_idx).RMSE        = RMSE_FED_over_time_tmp.'; % T×num_targets
end % End Monte-Carlo 

%%%%%%%%% RMSE struct (last sim)%%%%%%%%%
EKF_results_last                                    = monte_carlo_results(end).EKF_results;
federated_results_last                              = monte_carlo_results(end).federated_results;
RMSE_UAV_over_time                                  = zeros(num_uav, num_targets, T);
RMSE_FED_over_time                                  = zeros(num_targets, T);

for tgt_idx = 1:num_targets
    for t = 1:T
        gt_xy                                       = [true_traj(1,t,tgt_idx), true_traj(2,t,tgt_idx)];
        for uav_idx = 1:num_uav
            % est_xy                                  = EKF_results_last(uav_idx, tgt_idx).estimates(t,1:2);
            % RMSE_UAV_over_time(uav_idx, tgt_idx, t) = norm(gt_xy - est_xy);
            est_xyz = EKF_results_last(uav_idx, tgt_idx).estimates(t,1:3);
            RMSE_UAV_over_time(uav_idx, tgt_idx, t) = norm(gt_xyz - est_xyz);
        end
        % fed_xy                                      = federated_results_last(tgt_idx).estimates(t,1:2);
        % RMSE_FED_over_time(tgt_idx,t)               = norm(gt_xy - fed_xy);
        fed_xyz = federated_results_last(tgt_idx).estimates(t,1:3);
        RMSE_FED_over_time(tgt_idx,t) = norm(gt_xyz - fed_xyz);
    end
end

RMSE_UAV_overall                                    = squeeze(sqrt(mean(RMSE_UAV_over_time.^2, 3)) ); % [num_uav x num_targets]
RMSE_FED_overall                                    = sqrt(mean(RMSE_FED_over_time.^2, 2));            % [num_targets x 1]

rmse                                                = struct();
rmse.labels                                         = names;
if num_targets >= 1
    rmse.per_time.boat                              = RMSE_FED_over_time(1,:);   % 1×T
    rmse.overall.boat                               = RMSE_FED_overall(1);       % scalar
else
    rmse.per_time.boat                              = [];
    rmse.overall.boat                               = NaN;
end
if num_targets >= 2
    rmse.per_time.target                            = RMSE_FED_over_time(2,:); % 1×T
    rmse.overall.target                             = RMSE_FED_overall(2);     % scalar
else
    rmse.per_time.target                            = [];
    rmse.overall.target                             = NaN;
end
rmse.per_uav_overall                                = RMSE_UAV_overall; % [num_uav x num_targets]
rmse.time_vector                                    = dtv;

assignin('base','rmse',rmse);  % Output objective needed for DOE code

%% MC-averaged final-time RMSE per target (federated)
num_sims                                            = numel(monte_carlo_results);
T_final                                             = size(true_traj, 2);  % time index of final step
num_targets                                         = size(true_traj, 3);
RMSE_final_sq_accum                                 = zeros(num_targets,1);

for s = 1:num_sims
    if s <= numel(monte_carlo_RMSE_over_time) && ...
       isfield(monte_carlo_RMSE_over_time(s), 'RMSE') && ...
       ~isempty(monte_carlo_RMSE_over_time(s).RMSE)
        R_final = monte_carlo_RMSE_over_time(s).RMSE(end, :);   % 1×Ntgt
        RMSE_final_sq_accum = RMSE_final_sq_accum + (R_final.').^2;
    else
        fed_s = monte_carlo_results(s).federated_results;
        e_sq  = zeros(num_targets,1);
        for tgt_idx = 1:num_targets
            gt_xy   = [true_traj(1,T_final,tgt_idx), true_traj(2,T_final,tgt_idx)];
            est_xy  = fed_s(tgt_idx).estimates(T_final, 1:2);
            e_sq(tgt_idx) = sum((gt_xy - est_xy).^2);
        end
        RMSE_final_sq_accum = RMSE_final_sq_accum + e_sq;
    end
end

RMSE_final_MC = sqrt(RMSE_final_sq_accum / num_sims);   % [num_targets×1]

if exist('tr_crlb_AB', 'var') && exist('tr_crlb_AT', 'var')
    sqrt_tr_PCRLB = [sqrt(double(full(evalf(tr_crlb_AB)))); ...
        sqrt(double(full(evalf(tr_crlb_AT))))];
else
    sqrt_tr_PCRLB = [NaN; NaN];
end

fprintf('\n=== Monte-Carlo–Averaged Final-Time RMSE (St Dev) ===\n');
for tgt_idx = 1:num_targets
    if exist('names','var') && numel(names) >= tgt_idx
        fprintf('  %-16s : %8.3f m ->  sqrt(Tr(PCRLB)): %3.3f m\n', names{tgt_idx}, RMSE_final_MC(tgt_idx),sqrt_tr_PCRLB(tgt_idx));
    else
        fprintf('  Target %d          : %8.3f m\n', tgt_idx, RMSE_final_MC(tgt_idx));
    end
end

% output for DOE
assignin('base','RMSE_MC', RMSE_final_MC);

%% Plots
if data.plotfig && isequal(data.isDOE,false)
    plot_data = struct();
    plot_data.true_target_pos               = true_target_pos; % [XY,target,T]
    plot_data.num_uav                       = num_uav;
    plot_data.num_targets                   = num_targets;
    plot_data.num_time_steps                = T;
    plot_data.EKF_results                   = EKF_results_last;
    plot_data.federated_results             = federated_results_last;
    plot_data.initial_estimates             = initial_estimates;   % truth at t=1 (pre-noise)
    plot_data.tf                            = tf;
    plot_data.names                         = names;
    plot_data.monte_carlo_results           = monte_carlo_results;
    plot_data.monte_carlo_RMSE_over_time    = monte_carlo_RMSE_over_time;
    estimator_plots(plot_data);
end

%% Functions
function out = iff(cond, a, b)
if cond 
    out = a; 
else 
    out = b; 
end
end

function copyIfField(src, dstStruct, uavName, fld)
if isfield(src.(uavName), fld)
    dstStruct.(uavName).(fld) = src.(uavName).(fld);
end
end