function plotResults(data_nonCoop, results_nonCoop, data_Coop, results_Coop)
% PLOTRESULTS_CASADI  Plot trajectories and FOV using CasADi‑generated results

import casadi.*

% Convert non‑cooperative trajectories from MX → double
num_non = data_nonCoop.num_uav;
% for uavIdx = 1:num_non
%     idx = num2str(uavIdx);
%     xoutField = ['xOut' idx];
%     if isfield(results_nonCoop, xoutField) && ~isempty(results_nonCoop.(xoutField))
%         trajCas = extractTrajectories_nonCooperative( ...
%                       results_nonCoop.(xoutField), data_nonCoop);
%         fnames = fieldnames(trajCas.(['UAV' idx]));
%         for k = 1:numel(fnames)
%             f = fnames{k};
%             results_nonCoop.(['UAV' idx]).(f) = ...
%                 full(trajCas.(['UAV' idx]).(f));
%         end
%     end
% end

% Convert cooperative trajectories 
hasCoopData = (nargin == 4) && ~isempty(data_Coop) && ~isempty(results_Coop);
if hasCoopData && isfield(results_Coop, 'xout') && ~isempty(results_Coop.xout)
    trajCas = extractTrajectories_Cooperative( ...
                  results_Coop.xout, data_Coop);
    num_coop = data_Coop.num_uav;
    for uavIdx = 1:num_coop
        idx = num2str(uavIdx);
        fnames = fieldnames(trajCas.(['UAV' idx]));
        for k = 1:numel(fnames)
            f = fnames{k};
            results_Coop.(['UAV' idx]).(f) = ...
                full(trajCas.(['UAV' idx]).(f));
        end
    end
end

% Determine which dataset to plot
if hasCoopData
    data        = data_Coop;
    num_uav     = data.num_uav;
    BN          = data.BN;
    time        = data.time';
    alpha       = results_Coop.xout(end);
    tf          = data.tf * (1 + alpha);
    isGimballed = data.isGimballed;
else
    data        = data_nonCoop;
    num_uav     = data.num_uav;
    BN          = data.BN;
    time        = data.time';
    alpha       = results_nonCoop.xOut1(end);
    tf          = data.tf * (1 + alpha);
    isGimballed = data.isGimballed;
end

% Plot

% Plot NFZ & comm circles, boat & target
x_B   = data.X_B(1,:)';
y_B   = data.X_B(2,:)';
x_T   = data.x_T(1,:)';
y_T   = data.x_T(2,:)';
r_NFZ = data.r_NFZ;
rCom  = data.rCom;
theta = linspace(0,2*pi,360);

% UAV TRAJECTORIES (non-gimballed, num_uav > 1)
if ~isGimballed && (num_uav > 1)
    figure('Position',[300,300,1000,800]); hold on; grid on;
    % -- Non-Cooperative --
    for uavIdx = 1:num_uav
        idx = num2str(uavIdx);
        xout = ['xOut' idx];
        if ~isfield(results_nonCoop, xout) || isempty(results_nonCoop.(xout))
            fprintf('No non-coop trajectory for UAV %d; skipping.\n', uavIdx);
            continue;
        end
        % xA = BN * results_nonCoop.(['UAV' idx]).X_A;
        % yA = BN * results_nonCoop.(['UAV' idx]).Y_A;
        % plot(xA, yA, '--', 'DisplayName',['UAV ' idx ' Non‑coop']);
    end
    % -- Cooperative --
    try
        for uavIdx = 1:num_uav
            idx = num2str(uavIdx);
            X_A = BN * results_Coop.(['UAV' idx]).X_A;
            Y_A = BN * results_Coop.(['UAV' idx]).Y_A;
            plot(X_A, Y_A, '-', 'DisplayName',['UAV ' idx]);
        end

        fill(x_T(end)+r_NFZ*cos(theta), y_T(end)+r_NFZ*sin(theta), ...
             [1 0 0], 'FaceAlpha',0.1,'EdgeColor','r','LineStyle','--',...
             'DisplayName','No‑Fly Zone');
        fill(x_B(end)+rCom*cos(theta),  y_B(end)+rCom*sin(theta), ...
             [0.6 1 0.6], 'FaceAlpha',0.3,'EdgeColor','g','LineStyle','-',...
             'DisplayName','Comms Range');
        plot(x_B, y_B, 's','MarkerFaceColor','b','MarkerEdgeColor','b','DisplayName','Supporting USV');
        plot(x_T, y_T, 'd','MarkerFaceColor','r','MarkerEdgeColor','r','DisplayName','Target USV');
        legend('Location','bestoutside'); xlabel('x (m)'); ylabel('y (m)');
        axis equal; hold off;
    catch
        % ignore
    end

    % 1‑D Plots
    try
        figure('Position',[300,300,1000,900]);
        % create subplots
        for s = 1:4
            subplot(4,1,s); hold on; grid on;
            xlabel('Time (s)'); xlim([0 tf]);
        end
        ylabelTexts = {'x (m)','y (m)','\psi (deg)','\phi (deg)'};
        for s=1:4, subplot(4,1,s); ylabel(ylabelTexts{s}); end

        for uavIdx=1:num_uav
            idx = num2str(uavIdx);
            X_A   = BN * results_Coop.(['UAV' idx]).X_A;
            Y_A   = BN * results_Coop.(['UAV' idx]).Y_A;
            Psi_A = BN * results_Coop.(['UAV' idx]).Psi_A * (180/pi);
            Phi_A = BN * results_Coop.(['UAV' idx]).Phi_A * (180/pi);
            subplot(4,1,1), plot(time*tf, X_A,   'DisplayName',['UAV ' idx]);
            subplot(4,1,2), plot(time*tf, Y_A,   'DisplayName',['UAV ' idx]);
            subplot(4,1,3), plot(time*tf, Psi_A, 'DisplayName',['UAV ' idx]);
            subplot(4,1,4), plot(time*tf, Phi_A, 'DisplayName',['UAV ' idx]);
        end
        for s=1:4, subplot(4,1,s), legend('Location','bestoutside'); end
    catch
        % ignore
    end
end

% SINGLE FIXED‑FOV UAV (num_uav==1, not gimballed)
if ~isGimballed && (num_uav == 1)
    figure('Position',[300,300,1000,800]); hold on; grid on;
    % -- Non‑Coop skipped as before --
    % -- Coop Trajectory --
    try
        idx = '1';
        X_A = BN * results_Coop.UAV1.X_A;
        Y_A = BN * results_Coop.UAV1.Y_A;
        plot(X_A, Y_A, 'DisplayName','Fixed‑FOV UAV');
        % plot NFZ & comm, USV, target
        x_B = data.X_B(1,:); y_B = data.X_B(2,:);
        x_T = data.x_T(1,:); y_T = data.x_T(2,:);
        r_NFZ = data.r_NFZ; rCom = data.rCom;
        theta = linspace(0,2*pi,360);
        fill(x_T(end)+r_NFZ*cos(theta), y_T(end)+r_NFZ*sin(theta), ...
             [1 0 0], 'FaceAlpha',0.1,'EdgeColor','r','LineStyle','--',...
             'DisplayName','No‑Fly Zone');
        fill(x_B(end)+rCom*cos(theta), y_B(end)+rCom*sin(theta), ...
             [0.6 1 0.6], 'FaceAlpha',0.3,'EdgeColor','g','LineStyle','-',...
             'DisplayName','Comms Range');
        plot(x_B, y_B, 's','MarkerFaceColor','b','MarkerEdgeColor','b','DisplayName','Supporting USV');
        plot(x_T, y_T, 'd','MarkerFaceColor','r','MarkerEdgeColor','r','DisplayName','Target USV');
        legend('Location','bestoutside'); xlabel('x (m)'); ylabel('y (m)');
        axis equal; hold off;
    catch
        % ignore
    end

    % 1‑D Plots
    try
        figure('Position',[300,300,1000,800]);
        idx = '1';
        X_A   = BN * results_Coop.UAV1.X_A;
        Y_A   = BN * results_Coop.UAV1.Y_A;
        Psi_A = BN * results_Coop.UAV1.Psi_A*(180/pi);
        Phi_A = BN * results_Coop.UAV1.Phi_A*(180/pi);
        subplot(4,1,1); plot(time*tf, X_A,'DisplayName','x'); grid on; ylabel('x (m)'); xlim([0 tf]);
        subplot(4,1,2); plot(time*tf, Y_A,'DisplayName','y'); grid on; ylabel('y (m)'); xlim([0 tf]);
        subplot(4,1,3); plot(time*tf, Psi_A,'DisplayName','\psi'); grid on; ylabel('Yaw (deg)'); xlim([0 tf]);
        subplot(4,1,4); plot(time*tf, Phi_A,'DisplayName','\phi'); grid on; ylabel('Bank (deg)'); xlabel('Time (s)'); xlim([0 tf]);
    catch
        % ignore
    end

    % FOV Plot for Single Fixed‑FOV
    try
        idx = '1';
        X_A = BN * results_Coop.UAV1.X_A;
        Y_A = BN * results_Coop.UAV1.Y_A;
        Z_A = BN * results_Coop.UAV1.Z_A;
        camPos = [X_A.'; Y_A.'; Z_A.'];
        targPos = [data.x_T(1,end); data.x_T(2,end); 0];
        boatPos = [data.X_B(1,end); data.X_B(2,end); 0];
        traj = struct();  % build minimal traj struct for InFOV
        traj.UAV1.X_A       = X_A;
        traj.UAV1.Y_A       = Y_A;
        traj.UAV1.Psi_A     = BN*results_Coop.UAV1.Psi_A;
        traj.UAV1.tan_Phi_A = BN*results_Coop.UAV1.tan_Phi_A;
        traj.UAV1.Phi_A     = atan(traj.UAV1.tan_Phi_A);
        traj.UAV1.Z_A       = Z_A;
        inFOV_target  = InFOV(camPos, targPos, traj, idx, data);
        inFOV_boat    = InFOV(camPos, boatPos, traj, idx, data);
        resultsFOV.time = data.time(1:numel(X_A))*tf;
        resultsFOV.UAV_positions = [X_A, Y_A, Z_A];
        resultsFOV.target_in_view = [double(inFOV_target>.01).'; double(inFOV_boat>.01).'];
        targets = [data.x_T(1,end), data.x_T(2,end), 0;
                   data.X_B(1,end), data.X_B(2,end), 0];
        plotFOV(resultsFOV, targets, data, str2double(idx));
    catch ME
        fprintf('FOV plot failed: %s\n', ME.message);
    end
end


%%%%%%%%%%%%%%% Gimballed plots  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isGimballed && (num_uav == 1)
    %%%%%%%%%%%%%%%% Plot 2D trajectory
    try
        figure('Position',[300,300,1000,800]); hold on; grid on;
        idx = '1';

        X_A = BN * results_Coop.UAV1.X_A;
        Y_A = BN * results_Coop.UAV1.Y_A;
        plot(X_A, Y_A, 'DisplayName','Gimballed UAV');

        % Plot NFZ & comm, USV, target
        fill(x_T(end)+r_NFZ*cos(theta), y_T(end)+r_NFZ*sin(theta), ...
             [1 0 0], 'FaceAlpha',0.1,'EdgeColor','r','LineStyle','--',...
             'DisplayName','No‑Fly Zone');
        fill(x_B(end)+rCom*cos(theta), y_B(end)+rCom*sin(theta), ...
             [0.6 1 0.6], 'FaceAlpha',0.3,'EdgeColor','g','LineStyle','-',...
             'DisplayName','Comms Range');
        plot(x_B, y_B, 's','MarkerFaceColor','b','MarkerEdgeColor','b','DisplayName','Supporting USV');
        plot(x_T, y_T, 'd','MarkerFaceColor','r','MarkerEdgeColor','r','DisplayName','Target USV');
        legend('Location','bestoutside'); xlabel('x (m)'); ylabel('y (m)');
        axis equal; hold off;
    catch
        % ignore
    end

    %%%%%%%% Plot 1D state and gimbal angles
    try
        figure('Position',[300,300,1000,900]);
        for s = 1:6
            subplot(6,1,s); hold on; grid on;
            xlabel('Time (s)'); xlim([0 tf]);
        end
        ylabelTexts = {'x (m)','y (m)','\psi (deg)','\phi (deg)','\psi_g (deg)','\theta_g (deg)'};
        for s = 1:6
            subplot(6,1,s); ylabel(ylabelTexts{s});
        end
        
        X_A    = BN * results_Coop.UAV1.X_A;
        Y_A    = BN * results_Coop.UAV1.Y_A;
        Psi_A  = BN * results_Coop.UAV1.Psi_A * (180/pi);

        if isfield(results_Coop.UAV1, 'Phi_A') && ~isempty(results_Coop.UAV1.Phi_A)
            Phi_A = BN * results_Coop.UAV1.Phi_A * (180/pi);
        elseif isfield(results_Coop.UAV1, 'tan_Phi_A') && ~isempty(results_Coop.UAV1.tan_Phi_A)
            Phi_A = atan(BN * results_Coop.UAV1.tan_Phi_A) * (180/pi);
        else
            Phi_A = zeros(size(X_A));
        end

        cPsi_g   = BN * results_Coop.UAV1.cPsi_g_A;
        sPsi_g   = BN * results_Coop.UAV1.sPsi_g_A;
        cTheta_g = BN * results_Coop.UAV1.cTheta_g_A;
        sTheta_g = BN * results_Coop.UAV1.sTheta_g_A;

        psi_g   = atan2(sPsi_g,   cPsi_g)   * (180/pi);
        theta_g = atan2(sTheta_g, cTheta_g) * (180/pi);

        t_plot = time * tf;
        subplot(6,1,1); plot(t_plot, X_A, 'DisplayName','x');
        subplot(6,1,2); plot(t_plot, Y_A, 'DisplayName','y');
        subplot(6,1,3); plot(t_plot, Psi_A, 'DisplayName','\psi');
        subplot(6,1,4); plot(t_plot, Phi_A, 'DisplayName','\phi');
        subplot(6,1,5); plot(t_plot, psi_g, 'DisplayName','\psi_g');
        subplot(6,1,6); plot(t_plot, theta_g, 'DisplayName','\theta_g');
        
        for s = 1:6
            subplot(6,1,s); 
            legend('Location','bestoutside');
        end
    catch
        % ignore
    end

    %%%%%%%%%% FOV plot for single gimballed UAV
    % try
    %     idx = '1';
    %     % Use non‑cooperative results for gimballed FOV plot
    %     X_A = BN * results_nonCoop.UAV1.X_A;
    %     Y_A = BN * results_nonCoop.UAV1.Y_A;
    %     Z_A = BN * results_nonCoop.UAV1.Z_A;
    %     camPos = [X_A.'; Y_A.'; Z_A.'];
    %     targPos = [data.x_T(1,end); data.x_T(2,end); 0];
    %     boatPos = [data.X_B(1,end); data.X_B(2,end); 0];
    %     traj = struct();
    %     traj.UAV1.X_A        = X_A;
    %     traj.UAV1.Y_A        = Y_A;
    %     traj.UAV1.Psi_A      = BN * results_nonCoop.UAV1.Psi_A;
    %     if isfield(results_nonCoop.UAV1, 'tan_Phi_A') && ~isempty(results_nonCoop.UAV1.tan_Phi_A)
    %         traj.UAV1.tan_Phi_A   = BN * results_nonCoop.UAV1.tan_Phi_A;
    %         traj.UAV1.Phi_A       = atan(traj.UAV1.tan_Phi_A);
    %     elseif isfield(results_nonCoop.UAV1, 'Phi_A') && ~isempty(results_nonCoop.UAV1.Phi_A)
    %         traj.UAV1.Phi_A   = BN * results_nonCoop.UAV1.Phi_A;
    %         traj.UAV1.tan_Phi_A = tan(traj.UAV1.Phi_A);
    %     else
    %         traj.UAV1.Phi_A = zeros(size(X_A));
    %         traj.UAV1.tan_Phi_A = tan(traj.UAV1.Phi_A);
    %     end
    %     traj.UAV1.Z_A        = Z_A;
    %     traj.UAV1.cPsi_g_A   = BN * results_nonCoop.UAV1.cPsi_g_A;
    %     traj.UAV1.sPsi_g_A   = BN * results_nonCoop.UAV1.sPsi_g_A;
    %     traj.UAV1.cTheta_g_A = BN * results_nonCoop.UAV1.cTheta_g_A;
    %     traj.UAV1.sTheta_g_A = BN * results_nonCoop.UAV1.sTheta_g_A;
    %     % Compute visibility
    %     inFOV_target  = InFOV(camPos, targPos, traj, idx, data);
    %     inFOV_boat    = InFOV(camPos, boatPos, traj, idx, data);
    %     resultsFOV.time = data.time(1:numel(X_A)) * tf;
    %     resultsFOV.UAV_positions = [X_A, Y_A, Z_A];
    %     targVis = full(evalf(inFOV_target)) > 0.01;
    %     boatVis = full(evalf(inFOV_boat))   > 0.01;
    %     resultsFOV.target_in_view = [double(targVis).'; double(boatVis).'];
    %     targets = [data.x_T(1,end), data.x_T(2,end), 0; data.X_B(1,end), data.X_B(2,end), 0];
    %     plotFOV(resultsFOV, targets, data, str2double(idx));
    % catch ME
    %     fprintf('FOV plot failed: %s\n', ME.message);
    % end
    try
        idx = '1';
        X_A = BN * results_Coop.UAV1.X_A;
        Y_A = BN * results_Coop.UAV1.Y_A;
        Z_A = BN * results_Coop.UAV1.Z_A;

        camPos  = [X_A.'; Y_A.'; Z_A.'];
        targPos = [data.x_T(1,end); data.x_T(2,end); 0];
        boatPos = [data.X_B(1,end); data.X_B(2,end); 0];

        % Build minimal traj struct consumed by InFOV
        traj = struct();
        traj.UAV1.X_A   = X_A;
        traj.UAV1.Y_A   = Y_A;
        traj.UAV1.Psi_A = BN * results_Coop.UAV1.Psi_A;

        if isfield(results_Coop.UAV1, 'tan_Phi_A') && ~isempty(results_Coop.UAV1.tan_Phi_A)
            traj.UAV1.tan_Phi_A = BN * results_Coop.UAV1.tan_Phi_A;
            traj.UAV1.Phi_A     = atan(traj.UAV1.tan_Phi_A);
        elseif isfield(results_Coop.UAV1, 'Phi_A') && ~isempty(results_Coop.UAV1.Phi_A)
            traj.UAV1.Phi_A     = BN * results_Coop.UAV1.Phi_A;
            traj.UAV1.tan_Phi_A = tan(traj.UAV1.Phi_A);
        else
            traj.UAV1.Phi_A     = zeros(size(X_A));
            traj.UAV1.tan_Phi_A = tan(traj.UAV1.Phi_A);
        end

        traj.UAV1.Z_A        = Z_A;
        traj.UAV1.cPsi_g_A   = BN * results_Coop.UAV1.cPsi_g_A;
        traj.UAV1.sPsi_g_A   = BN * results_Coop.UAV1.sPsi_g_A;
        traj.UAV1.cTheta_g_A = BN * results_Coop.UAV1.cTheta_g_A;
        traj.UAV1.sTheta_g_A = BN * results_Coop.UAV1.sTheta_g_A;

        % Visibility
        inFOV_target = InFOV(camPos, targPos, traj, idx, data);
        inFOV_boat   = InFOV(camPos, boatPos, traj, idx, data);

        resultsFOV.time          = data.time(1:numel(X_A)) * tf;
        resultsFOV.UAV_positions = [X_A, Y_A, Z_A];

        targVis = full((inFOV_target)) > 0.0001;
        boatVis = full((inFOV_boat))   > 0.0001;
        resultsFOV.target_in_view = [double(targVis).'; double(boatVis).'];

        targets = [data.x_T(1,end), data.x_T(2,end), 0; ...
                   data.X_B(1,end), data.X_B(2,end), 0];

        plotFOV(resultsFOV, targets, data, str2double(idx));
    catch ME
        fprintf('FOV plot failed: %s\n', ME.message);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotFOV(resultsFOV, targets, data, uavIdx)
% plotFOV  Color-codes UAV flight path according to targets in FOV, using enhanced color scheme.

% Extract positions and visibility matrix
UAVpos    = resultsFOV.UAV_positions;
visMat    = resultsFOV.target_in_view;
num_steps = size(UAVpos,1);
num_targets = size(targets,1);

% Define better, clearer colors
noColor    = [0.6 0.6 0.6];            % Light gray
targetColor = [1 0 0];                 % Bright red (Target USV)
cooperativeColor = [0 0 1];    % Bright blue (Supporting USV)
bothColor  = [0.494 0.184 0.556];       % Purple

% Colors array for quick reference
baseColors = [targetColor; cooperativeColor];

% Initialize color array
posColors = zeros(num_steps,3);

% Assign color to each time step
for k = 1:num_steps
    visible = find(visMat(:,k) == 1);
    if isempty(visible)
        posColors(k,:) = noColor;
    elseif length(visible) == 1
        posColors(k,:) = baseColors(visible,:);
    else
        posColors(k,:) = bothColor;
    end
end

% Create figure
if data.isGimballed
    figName = sprintf('Ideal Gimbal FOV (UAV %d)', uavIdx);
else
    figName = sprintf('Ideal Fixed-FOV (UAV %d)', uavIdx);
end

figure('Name', figName, 'NumberTitle','off','Position',[100 100 900 700]);
hold on; grid on;
scatter(UAVpos(:,1), UAVpos(:,2), 36, posColors, 'filled');

% Plot targets
tgtNames = {'Target USV','Supporting USV'};
legendHandles = gobjects(num_targets,1);

for iT = 1:num_targets
    if iT == 1
        marker = 'd'; edgeColor = 'r'; faceColor = 'r';
    else
        marker = 's'; edgeColor = 'b'; faceColor = 'b';
    end
    legendHandles(iT) = plot(targets(iT,1), targets(iT,2), marker, ...
        'MarkerEdgeColor',edgeColor,'MarkerFaceColor',faceColor,'LineWidth',2, ...
        'DisplayName', tgtNames{iT});
end

% Plot No-Fly Zone and Communication Range
theta = linspace(0,2*pi,100);
r_NFZ = data.r_NFZ;
rCom  = data.rCom;

xC_NFZ = targets(1,1) + r_NFZ*cos(theta);
yC_NFZ = targets(1,2) + r_NFZ*sin(theta);
L_NFZ = fill(xC_NFZ,yC_NFZ,[1 0 0], 'FaceAlpha', 0.1, ...
    'EdgeColor', 'r', 'LineStyle', '--', 'DisplayName', 'No-fly Zone');

xC_Com = targets(2,1) + rCom*cos(theta);
yC_Com = targets(2,2) + rCom*sin(theta);
L_com = fill(xC_Com,yC_Com,[0.6 1 0.6], 'FaceAlpha', 0.3, ...
    'EdgeColor', 'g', 'LineStyle', '-', 'DisplayName', 'Communication Range');

% Dummy points for FOV legend
dummyHandles = gobjects(4,1);
legendFOV = {'None in FOV', 'Target USV in FOV', 'Supporting USV in FOV', 'Both in FOV'};
dummyColors = [noColor; targetColor; cooperativeColor; bothColor];

for idx = 1:4
    dummyHandles(idx) = scatter(nan, nan, 36, dummyColors(idx,:), 'filled');
end

% Finalize legend
legend([legendHandles; L_NFZ; L_com; dummyHandles], ...
    [tgtNames(:); {'No-fly Zone'; 'Communication Range'}; legendFOV(:)], ...
    'Location','bestoutside');

xlabel('X Position (m)');
ylabel('Y Position (m)');
% title(figName);
axis equal;
hold off;
end

