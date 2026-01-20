% estimator_plots.m
% Produces figures for trajectories, covariance ellipses, RMSE over time, and histograms.
function estimator_plots(data)
true_target_pos   = data.true_target_pos;   % [2 x num_targets x T]
num_UAVs          = data.num_uav;
num_targets       = data.num_targets;
num_time_steps    = data.num_time_steps;
EKF_results       = data.EKF_results;           % single run (for plots 1–3)
initial_estimates = data.initial_estimates;
federated_results = data.federated_results;     % single run (for plots 1–3)
tf                = data.tf;
names             = data.names;

hasMC = isfield(data,'monte_carlo_results') && ~isempty(data.monte_carlo_results);

%% %%%%%% Position estimates (single simulation visuals) %%%%%%%%%%%
try
    figure('Name','EKF Position Estimates','Position',[300,300,1200,600]);
    colors = lines(num_UAVs);
    for target_idx = 1:num_targets
        subplot(1, num_targets, target_idx); hold on;
        for uav_idx = 1:num_UAVs
            X = EKF_results(uav_idx, target_idx).estimates(:,1);
            Y = EKF_results(uav_idx, target_idx).estimates(:,2);
            ct = createcolortable(uav_idx, num_UAVs, num_time_steps);
            plot(X(end), Y(end), 'o', 'Color', colors(uav_idx,:), 'LineWidth', 1, ...
                'DisplayName', sprintf('UAV %d',uav_idx));
            for k = 1:num_time_steps-1
                plot(X(k:k+1), Y(k:k+1), 'Color', ct(k,:), 'HandleVisibility','off','LineWidth',2);
            end
        end
        scatter(true_target_pos(1,target_idx,end), true_target_pos(2,target_idx,end), ...
            50,'r','diamond','LineWidth',2,'DisplayName','True Pos');
        scatter(initial_estimates(target_idx,1), initial_estimates(target_idx,2), ...
            80,'r','*','LineWidth',1.5,'DisplayName','Init (truth t=1)');
        legend(Location="best");
        xlabel('X (m)'); ylabel('Y (m)');
        title(sprintf('%s Localization', names{target_idx}));
        grid on; axis equal; hold off;
    end
catch ME
    warning('PositionEstimates: %s', ME.message); %#ok<MEXCEP>
end

% %% Individual covariance ellipses for target
% try
%     if num_targets >= 2
%         target_idx = 2;
%         figure('Name',sprintf('Individual UAV Covariance Ellipses for %s',names{target_idx}),'Position',[300,300,1000,800]); hold on;
% 
%         for uav_idx = 1:num_UAVs
%             plot_uncertainty_ellipse(EKF_results(uav_idx,target_idx).estimates(end,1:2), ...
%                 EKF_results(uav_idx,target_idx).covariances(1:2,1:2,end), uav_idx);
%         end
%         scatter(true_target_pos(1,target_idx,end), true_target_pos(2,target_idx,end), ...
%             50,'r','diamond','LineWidth',2,'DisplayName','True Pos');
%         scatter(initial_estimates(target_idx,1), initial_estimates(target_idx,2), ...
%             80,'r','*','LineWidth',1.5,'DisplayName','Init');
%         legend(Location="best"); xlabel('X (m)'); ylabel('Y (m)');
%         title(sprintf('Individual UAV Covariance Ellipses for %s at Final Time',names{target_idx})); grid on; axis equal; hold off;
%     end
% catch ME
%     warning('IndividualEllipses: %s', ME.message);%#ok<MEXCEP>
% end

%% Federated covariance ellipse for target
try
    if num_targets >= 2
        target_idx = 2;

        %==== Main figure/axes ====
        figH = figure('Name',sprintf('Federated Covariance Ellipse for %s',names{target_idx}), ...
            'Position',[300,300,1000,800]);
        ax = axes(figH); hold(ax, 'on');

        est_avg = federated_results(target_idx).estimates(end,1:2);
        plot_uncertainty_ellipse( est_avg, federated_results(target_idx).covariances(1:2,1:2,end), 0);
        % scatter(ax, est_avg(1), est_avg(2), 36, 'k', 'filled', ...
        %         'DisplayName','Federated Center');
        scatter(ax, true_target_pos(1,target_idx,end), true_target_pos(2,target_idx,end), ...
            50,'r','diamond','LineWidth',2,'DisplayName','True Pos');

        scatter(ax, initial_estimates(target_idx,1), initial_estimates(target_idx,2), ...
            80,'r','*','LineWidth',1.5,'DisplayName','Init');

        legend(ax, "Location","best");
        xlabel(ax, 'X (m)'); ylabel(ax, 'Y (m)');
        title(ax, sprintf('Federated Covariance Ellipse for %s at Final Time',names{target_idx}));
        grid(ax, 'on'); axis(ax, 'equal');
        hold(ax, 'off');

        % Inset zoom
        inset_size = struct('w',0.25, 'h',0.25);  % normalized (figure) units
        margin     = 0.02;                        % normalized
        weight_cfg = struct('use_weights', true, 'special_w', 5.0);

        % Protected points
        prot_pts = [ ...
            est_avg(:)'; ...
            true_target_pos(1,target_idx,end), true_target_pos(2,target_idx,end); ...
            initial_estimates(target_idx,1),   initial_estimates(target_idx,2) ...
            ];

        best_pos = local_best_inset_position(ax, inset_size, margin, prot_pts, weight_cfg);

        inset_ax = axes('Parent', figH, 'Position', [best_pos.x, best_pos.y, inset_size.w, inset_size.h]);
        box(inset_ax, 'on'); hold(inset_ax, 'on');

        P = federated_results(target_idx).covariances(1:2,1:2,end);
        plot_uncertainty_ellipsoid(est_avg, P, 0, 'k', inset_ax);
        scatter(inset_ax, est_avg(1), est_avg(2), 'ko','filled');
        scatter(inset_ax, true_target_pos(1,target_idx,end), true_target_pos(2,target_idx,end), ...
            25,'r','diamond','LineWidth',1.2);

        [~, D] = eig(P);
        stds = sqrt(diag(D));
        zf = 1.5;
        xlim(inset_ax, est_avg(1) + zf*[-stds(1), stds(1)]);
        ylim(inset_ax, est_avg(2) + zf*[-stds(2), stds(2)]);
        axis(inset_ax, 'equal'); grid(inset_ax, 'on'); hold(inset_ax, 'off');
    end

catch ME
    warning('FederatedEllipse: %s', ME.message); %#ok<MEXCEP>
end


%% RMSE over time (Monte-Carlo average)
try
    % Time vector for plotting
    tt = tf/num_time_steps*(1:num_time_steps);

    % Check for precomputed per-time federated RMSE curves
    hasMCcurves = hasMC && isfield(data,'monte_carlo_RMSE_over_time') && ...
        ~isempty(data.monte_carlo_RMSE_over_time) && ...
        isfield(data.monte_carlo_RMSE_over_time(1), 'RMSE');

    if hasMC
        num_sims = numel(data.monte_carlo_results);
    else
        num_sims = 0;
    end

    for target_idx = 1:num_targets

        % ----- Federated, MC-averaged per-time RMSE -----
        if hasMCcurves
            acc_fed_sq = zeros(1, num_time_steps);
            for s = 1:num_sims
                R = data.monte_carlo_RMSE_over_time(s).RMSE; % T×Ntgt
                r_target = R(:,target_idx).';                % 1×T
                acc_fed_sq = acc_fed_sq + r_target.^2;
            end
            RMSE_FED_MC = sqrt(acc_fed_sq / num_sims);       % 1×T
        elseif hasMC
            acc_fed_sq = zeros(1, num_time_steps);
            for s = 1:num_sims
                fed_s = data.monte_carlo_results(s).federated_results;
                r2    = zeros(1, num_time_steps);
                for t = 1:num_time_steps
                    gt_xy  = [true_target_pos(1,target_idx,t), true_target_pos(2,target_idx,t)];
                    est_xy = fed_s(target_idx).estimates(t,1:2);
                    r2(t)  = sum((gt_xy - est_xy).^2);
                end
                acc_fed_sq = acc_fed_sq + r2;
            end
            RMSE_FED_MC = sqrt(acc_fed_sq / num_sims);       % 1×T
        else
            RMSE_FED_MC = zeros(1, num_time_steps);
            for t = 1:num_time_steps
                gt_xy  = [true_target_pos(1,target_idx,t), true_target_pos(2,target_idx,t)];
                est_xy = federated_results(target_idx).estimates(t,1:2);
                RMSE_FED_MC(t) = norm(gt_xy - est_xy);
            end
        end

 % Per-UAV curves
            RMSE_UAV_plot = zeros(num_UAVs, num_time_steps);
            if hasMC
                for uav_idx = 1:num_UAVs
                    acc_sq = zeros(1, num_time_steps);
                    for s = 1:num_sims
                        EKF_s = data.monte_carlo_results(s).EKF_results;
                        r2    = zeros(1, num_time_steps);
                        for t = 1:num_time_steps
                            gt_xy  = [true_target_pos(1,target_idx,t), true_target_pos(2,target_idx,t)];
                            est_xy = EKF_s(uav_idx, target_idx).estimates(t,1:2);
                            r2(t)  = sum((gt_xy - est_xy).^2);
                        end
                        acc_sq = acc_sq + r2;
                    end
                    RMSE_UAV_plot(uav_idx,:) = sqrt(acc_sq / num_sims);
                end
            else
                for uav_idx = 1:num_UAVs
                    for t = 1:num_time_steps
                        gt_xy  = [true_target_pos(1,target_idx,t), true_target_pos(2,target_idx,t)];
                        est_xy = EKF_results(uav_idx, target_idx).estimates(t,1:2);
                        RMSE_UAV_plot(uav_idx,t) = norm(gt_xy - est_xy);
                    end
                end
            end

            figure('Name',sprintf('RMSE Over Time – %s', names{target_idx}), ...
                   'Position',[300,300,1000,600]); hold on;
            for uav_idx = 1:num_UAVs
                plot(tt, RMSE_UAV_plot(uav_idx,:), 'LineWidth', 1.5, ...
                    'DisplayName', sprintf('UAV %d', uav_idx));
            end
            legend('Location','best'); xlabel('Time (s)'); ylabel('RMSE (m)');
            title(sprintf('%s RMSE%s', names{target_idx}, ternary(hasMC,' (Monte-Carlo Average)','')));
            grid on; hold off;
        end
    catch ME
        warning('RMSEPlot: %s', ME.message); %#ok<MEXCEP>
    end

%% MC-Averaged Per-UAV Filter Uncertainty over Time 
try
        if ~hasMC
            warning('UncertaintyPlot: Monte-Carlo containers not found; skipping MC-averaged uncertainty.');
        else
            num_sims = numel(data.monte_carlo_results);
            tt = tf/num_time_steps*(1:num_time_steps);

            for target_idx = 1:num_targets
                sigma_r_all = zeros(num_UAVs, num_time_steps);
                for uav_idx = 1:num_UAVs
                    for t = 1:num_time_steps
                        accP = zeros(2,2);
                        for s = 1:num_sims
                            Pfull = data.monte_carlo_results(s) ...
                                        .EKF_results(uav_idx, target_idx) ...
                                        .covariances(:,:,t);
                            accP = accP + Pfull(1:2,1:2);
                        end
                        Pbar = accP / num_sims;
                        trP  = max(trace((Pbar+Pbar.')/2), 0);
                        sigma_r_all(uav_idx,t) = sqrt(trP);
                    end
                end

                figure('Name', sprintf('MC-Averaged Filter Uncertainty – %s', names{target_idx}), ...
                       'Position', [320,320,1000,600]); hold on;
                for uav_idx = 1:num_UAVs
                    plot(tt, sigma_r_all(uav_idx,:), 'LineWidth', 1.6, ...
                         'DisplayName', sprintf('UAV %d', uav_idx));
                end
                xlabel('Time (s)'); ylabel('Uncertainty (m)'); grid on;
                title(sprintf('%s Filter Uncertainty (MC Average, radial 1\\sigma)', names{target_idx}));
                legend('Location','best'); hold off;
            end
        end
    catch ME
        warning('UncertaintyPlot: %s', ME.message); %#ok<MEXCEP>
    end

% Overlaid Final-Time RMSE Histogram (per UAV + federated)
    %    (Requires Monte-Carlo runs)
    try
        if ~hasMC
            warning('Histogram: Monte-Carlo containers not found; skipping final-time histogram.');
            return;
        end

        num_sims = numel(data.monte_carlo_results);
        colors_all = lines(num_UAVs + 1);   % +1 for federated

        % Choose which target to show (mirrors your other figures)
        tgt_list = 1:num_targets;
        if num_targets >= 2
            tgt_list = 2;    % show Target USV only (common in your figures)
        end

        for target_idx = tgt_list
            % Collect final-time RMSEs over sims
            T = size(true_target_pos, 3);
            final_err = zeros(num_sims, num_UAVs + 1); % [sims x (UAVs + federated)]

            for s = 1:num_sims
                EKF_s = data.monte_carlo_results(s).EKF_results;
                FED_s = data.monte_carlo_results(s).federated_results;

                gt_xy = [ true_target_pos(1,target_idx,T), true_target_pos(2,target_idx,T) ];

                % Each UAV
                for uav_idx = 1:num_UAVs
                    est_xy = EKF_s(uav_idx, target_idx).estimates(T,1:2);
                    final_err(s, uav_idx) = norm(gt_xy - est_xy);
                end
                % Federated
                est_xy_fed = FED_s(target_idx).estimates(T,1:2);
                final_err(s, num_UAVs+1) = norm(gt_xy - est_xy_fed);
            end

            % Plot overlaid histograms
            figure('Name', sprintf('Final-Time RMSE Histogram – %s', names{target_idx}), ...
                   'Position', [100,100,1000,800]);
            hold on; grid on;

            % Common bins
            bin_width = 0.125;                         % meters
            min_val   = min(final_err(:));
            max_val   = max(final_err(:));
            if min_val == max_val
                % degenerate case (all the same) – widen slightly
                min_val = min_val - 0.5; max_val = max_val + 0.5;
            end
            bin_edges = min_val:bin_width:max_val;

            % Names for legend
            path_names = [ arrayfun(@(i) sprintf('UAV %d',i), 1:num_UAVs, 'uni',0), {'Federated'} ];

            % Draw histograms
            for k = 1:(num_UAVs+1)
                histogram(final_err(:,k), bin_edges, ...
                    'FaceColor', colors_all(k,:), ...
                    'FaceAlpha', 0.35, ...
                    'EdgeColor', 'none', ...
                    'DisplayName', sprintf('%s RMSE', path_names{k}));
            end

            % Mean lines
            for k = 1:(num_UAVs+1)
                mval = mean(final_err(:,k));
                xline(mval, '--', 'LineWidth', 1.5, 'Color', colors_all(k,:), ...
                    'DisplayName', sprintf('%s mean', path_names{k}));
            end

            xlabel('Final-time RMSE (m)'); ylabel('Monte-Carlo Trials');
            title(sprintf('Overlaid Final-Time RMSE Histogram – %s', names{target_idx}));
            legend('Location','northeast');
            hold off;
        end
    catch ME
        warning('Histogram: %s', ME.message);
    end  

% 6) MC-Averaged Final Location & Averaged Covariance Ellipses
try
    if ~hasMC
        warning('MCavgEllipse: Monte-Carlo results not found; using single run instead.');
    end

    % total number of MC sims (or pretend it's 1 if none)
    if hasMC
        num_sims = numel(data.monte_carlo_results);
    else
        num_sims = 1;
    end

    % final time index (same T used elsewhere in this file)
    Tfinal = num_time_steps;

    for target_idx = 1:num_targets
        % Prepare figure for this target
        figName = sprintf('MC-Averaged Final State & Avg Covariance – %s', names{target_idx});
        figure('Name', figName, 'Position', [300, 300, 1000, 800]); hold on; grid on;
        title(figName); xlabel('X (m)'); ylabel('Y (m)'); axis equal;

        % ===== Federated MC average =====
        if hasMC
            mu_stack = zeros(2, num_sims);
            P_sum    = zeros(2,2);
            for s = 1:num_sims
                fed_s = data.monte_carlo_results(s).federated_results;
                mu    = fed_s(target_idx).estimates(Tfinal,1:2).';        % 2x1
                Pxy   = fed_s(target_idx).covariances(1:2,1:2,Tfinal);    % 2x2
                mu_stack(:,s) = mu;
                P_sum = P_sum + (Pxy + Pxy.' )/2;
            end
            mu_bar_fed = mean(mu_stack, 2);       % 2x1
            P_bar_fed  = P_sum / num_sims;        % 2x2  (simple average of covariances)
        else
            % single-run fallback
            mu_bar_fed = data.federated_results(target_idx).estimates(Tfinal,1:2).';
            P_bar_fed  = data.federated_results(target_idx).covariances(1:2,1:2,Tfinal);
        end

        % Draw federated average ellipse (uav_idx=0 => black)
        plot_uncertainty_ellipse(mu_bar_fed.', P_bar_fed, 0);

        % True & init markers (for context)
        scatter(true_target_pos(1,target_idx,end), true_target_pos(2,target_idx,end), ...
            50,'r','diamond','LineWidth',2,'DisplayName','True Pos');
        scatter(initial_estimates(target_idx,1), initial_estimates(target_idx,2), ...
            80,'r','*','LineWidth',1.5,'DisplayName','Init');

        % ===== Per-UAV MC averages =====
        for uav_idx = 1:num_UAVs
            if hasMC
                mu_stack = zeros(2, num_sims);
                P_sum    = zeros(2,2);
                for s = 1:num_sims
                    EKF_s = data.monte_carlo_results(s).EKF_results;
                    mu    = EKF_s(uav_idx, target_idx).estimates(Tfinal,1:2).';       % 2x1
                    Pxy   = EKF_s(uav_idx, target_idx).covariances(1:2,1:2,Tfinal);   % 2x2
                    mu_stack(:,s) = mu;
                    P_sum = P_sum + (Pxy + Pxy.')/2;
                end
                mu_bar = mean(mu_stack, 2);
                P_bar  = P_sum / num_sims;
            else
                % single-run fallback
                mu_bar = data.EKF_results(uav_idx, target_idx).estimates(Tfinal,1:2).';
                P_bar  = data.EKF_results(uav_idx, target_idx).covariances(1:2,1:2,Tfinal);
            end

            % Draw the averaged ellipse for this UAV (colored by index)
            plot_uncertainty_ellipse(mu_bar.', P_bar, uav_idx);
        end

        legend('Location','best');
        hold off;
    end
catch ME
    warning('MCavgEllipse: %s', ME.message);
end


end


%% %%%%% Functions %%%%%%
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function ct = createcolortable(uav_idx, num_uav, t)
    col = lines(num_uav);
    colorA = 2.*col(uav_idx,:); colorA(colorA>1) = 1;
    colorB = 0.25.*colorA;
    ct = interp1([1/t 1],[colorA; colorB],(1:t)/t);
end


function plot_uncertainty_ellipse(position, covariance, uav_idx)
    k     = 1; % number of st.dev. to plot
    th    = linspace(0, 2*pi, 100);
    circ  = [cos(th); sin(th)];
    [V,D] = eig(covariance);
    E     = V*(k*sqrt(D))*circ;
    E     = E + position(:);
    if uav_idx == 0
        plot(E(1,:), E(2,:), 'k', 'LineWidth', 2, 'DisplayName', 'Federated');
    else
        colorMap = lines(10);
        plot(E(1,:), E(2,:), 'Color', colorMap(uav_idx,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('UAV %d', uav_idx));
        plot(position(1), position(2), '.', 'MarkerSize', 12, 'HandleVisibility','off');
    end
end


function plot_uncertainty_ellipsoid(position, covariance, uav_idx, colorval, ax)
    if nargin < 5, ax = gca; end
    [V, D] = eig(covariance);
    th     = linspace(0, 2*pi, 100);
    E      = V*sqrt(D)*[cos(th); sin(th)];
    E      = bsxfun(@plus, E, position(:));
    if isnumeric(colorval), col = colorval; else, col = [0 0 0]; end
    if uav_idx == 0
        plot(ax, E(1,:), E(2,:), 'LineWidth', 1, 'Color', col, 'DisplayName','Federated');
        plot(ax, position(1), position(2),'.','MarkerSize',8,'LineWidth',2,'Color',col,'HandleVisibility','off');
    else
        plot(ax, E(1,:), E(2,:), 'LineWidth', 1, 'Color', col, 'DisplayName', sprintf('UAV %d', uav_idx));
        plot(ax, position(1), position(2),'.','MarkerSize',8,'LineWidth',2,'Color',col,'HandleVisibility','off');
    end
end


function best_pos = local_best_inset_position(ax, inset_size, margin, prot_pts, weight_cfg)
    arguments
        ax (1,1) matlab.graphics.axis.Axes
        inset_size struct
        margin (1,1) double {mustBeNonnegative} = 0.01
        prot_pts double = zeros(0,2)
        weight_cfg struct = struct('use_weights',false,'special_w',1.0)
    end

    fig = ancestor(ax,'figure');
    oldUnitsFig = fig.Units; oldUnitsAx = ax.Units;
    fig.Units = 'normalized'; ax.Units = 'normalized';

    axpos = ax.Position; % [x y w h]
    xlimv = xlim(ax); ylimv = ylim(ax);
    dx    = xlimv(2) - xlimv(1);
    dy    = ylimv(2) - ylimv(1);

    frac_w = inset_size.w / axpos(3);
    frac_h = inset_size.h / axpos(4);

    cover_dx = frac_w * dx;
    cover_dy = frac_h * dy;

    cand(1).name = 'TR'; cand(1).x = axpos(1) + axpos(3) - inset_size.w - margin; cand(1).y = axpos(2) + axpos(4) - inset_size.h - margin;
    cand(2).name = 'TL'; cand(2).x = axpos(1) + margin;                                cand(2).y = axpos(2) + axpos(4) - inset_size.h - margin;
    cand(3).name = 'BR'; cand(3).x = axpos(1) + axpos(3) - inset_size.w - margin;      cand(3).y = axpos(2) + margin;
    cand(4).name = 'BL'; cand(4).x = axpos(1) + margin;                                cand(4).y = axpos(2) + margin;

    for k = 1:numel(cand)
        cand(k).x = max(axpos(1)+margin, min(cand(k).x, axpos(1)+axpos(3)-inset_size.w-margin));
        cand(k).y = max(axpos(2)+margin, min(cand(k).y, axpos(2)+axpos(4)-inset_size.h-margin));
    end

    ch = ax.Children;
    Xall = []; Yall = [];
    for i = 1:numel(ch)
        if isprop(ch(i),'XData') && isprop(ch(i),'YData') && strcmpi(ch(i).Visible,'on')
            Xi = ch(i).XData; Yi = ch(i).YData;
            try 
                Xi = Xi(:); Yi = Yi(:); 
            catch 
                Xi = []; Yi = []; 
            end
            v = ~isnan(Xi) & ~isnan(Yi);
            Xall = [Xall; Xi(v)]; %#ok<AGROW>
            Yall = [Yall; Yi(v)]; %#ok<AGROW>
        end
    end

    costs = zeros(1, numel(cand));
    for k = 1:numel(cand)
        switch cand(k).name
            case 'TR', xR = [xlimv(2)-cover_dx, xlimv(2)]; yR = [ylimv(2)-cover_dy, ylimv(2)];
            case 'TL', xR = [xlimv(1), xlimv(1)+cover_dx]; yR = [ylimv(2)-cover_dy, ylimv(2)];
            case 'BR', xR = [xlimv(2)-cover_dx, xlimv(2)]; yR = [ylimv(1), ylimv(1)+cover_dy];
            case 'BL', xR = [xlimv(1), xlimv(1)+cover_dx]; yR = [ylimv(1), ylimv(1)+cover_dy];
        end
        in_main = (Xall >= xR(1)) & (Xall <= xR(2)) & (Yall >= yR(1)) & (Yall <= yR(2));
        costs(k) = nnz(in_main);
        if weight_cfg.use_weights && ~isempty(prot_pts)
            px = prot_pts(:,1); py = prot_pts(:,2);
            in_prot = (px >= xR(1)) & (px <= xR(2)) & (py >= yR(1)) & (py <= yR(2));
            costs(k) = costs(k) + weight_cfg.special_w * nnz(in_prot);
        end
    end
    [~, idx] = min(costs(:));
    best_pos = struct('x', cand(idx).x, 'y', cand(idx).y);

    fig.Units = oldUnitsFig; ax.Units = oldUnitsAx;
end