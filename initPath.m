function [Cx, Cu] = initPath(data, uavIdx)
% INITPATH_FEASIBLE Generates a kinematically feasible initial guess.
% Uses Cubic Bezier curves for smooth kinematic connections and fits
% the Mixed Bernstein-Fourier basis to the result.
%
% FIXED: Dimension mismatch in evaluateBezierPhysics (enforced row vectors).

% -------------------------------------------------------------------------
% 1. Setup & Unpacking
% -------------------------------------------------------------------------
idx_str = num2str(uavIdx);
x0      = data.(['x0_A' idx_str]); % [x, y, psi]
xf      = data.(['xf_A' idx_str]); % [x, y, psi]
V_A     = data.V_A;
tf      = data.tf;
Nx      = data.Nx;
Nu      = data.Nu;
N_tot   = data.N_tot;
g       = data.g;
phi_lim = data.control_limit;

% Define dense grid for physics generation
n_fit = 400; 
t_fit = linspace(0, 1, n_fit)';

% -------------------------------------------------------------------------
% 2. Solve Geometrically Feasible Path
% -------------------------------------------------------------------------
% We need a path of length L_target ~ V_A * tf to satisfy time constraints.
L_target  = V_A * tf;

% FORCE ROW VECTORS [1x2]
p_start   = reshape(x0(1:2), 1, 2);
p_end     = reshape(xf(1:2), 1, 2);

psi_start = x0(3);

% Determine End Heading
if numel(xf) >= 3
    psi_end = xf(3);
else
    % If no end heading specified, point smoothly towards the target
    delta = p_end - p_start;
    psi_end = atan2(delta(2), delta(1)); 
end

% Bisection Search for Bezier Control Point Distance (R)
R_min  = 0; 
R_max  = L_target * 5; 
R_best = norm(p_end - p_start)/2;

for iter = 1:20
    R_curr = (R_min + R_max) / 2;
    L_curr = computeBezierLength(p_start, psi_start, p_end, psi_end, R_curr);
    
    if abs(L_curr - L_target) < 2.0 % Tolerance (meters)
        R_best = R_curr;
        break;
    end
    
    if L_curr < L_target
        R_min = R_curr; 
    else
        R_max = R_curr; 
    end
    R_best = R_curr;
end

% -------------------------------------------------------------------------
% 3. Generate Physics (Pos, Vel, Psi, Phi)
% -------------------------------------------------------------------------
[pos, ~, psi, phi] = evaluateBezierPhysics(p_start, psi_start, p_end, psi_end, ...
                                           R_best, t_fit, V_A, g, phi_lim);

% -------------------------------------------------------------------------
% 4. Fit Basis Coefficients
% -------------------------------------------------------------------------
% Construct Local Basis (Must match Main Script definitions exactly)
BN  = bernsteinMatrix_local(data.N_B, t_fit);
FN  = FourierMatrix_local(data.N_F, t_fit);
FBN = [BN, FN]; % [n_fit x N_tot]

% Regularized Least Squares
reg = 1e-10; 
Cx  = zeros(N_tot, Nx);
Cu  = zeros(N_tot, Nu);

% A. Fit Position [x, y]
Cx(:,1) = ridgeLS(FBN, pos(:,1), reg);
Cx(:,2) = ridgeLS(FBN, pos(:,2), reg);

% B. Fit Heading [psi] (Unwrap to avoid jumps)
Cx(:,3) = ridgeLS(FBN, unwrap(psi), reg);

% C. Fit Gimbal States (if applicable)
if Nx > 3
    % Simple heuristic: Point Gimbal at Target 1 constantly
    xT_cols = data.x_T; % [3 x N]
    
    % Interpolate Target Trajectory to match fit grid
    % Assuming x_T is defined on data.time_opt (0 to 1)
    xT_dense = zeros(3, n_fit);
    for k = 1:3
        xT_dense(k,:) = interp1(linspace(0,1,size(xT_cols,2)), xT_cols(k,:), t_fit);
    end
    
    % Vectors
    dx = xT_dense(1,:)' - pos(:,1);
    dy = xT_dense(2,:)' - pos(:,2);
    h  = data.(['z_A' idx_str]);
    dz = xT_dense(3,:)' - h;
    
    % Angles
    az_global = atan2(dy, dx);
    psi_g     = az_global - psi; % Relative Pan
    
    dist_2d   = sqrt(dx.^2 + dy.^2);
    theta_g   = atan2(dz, dist_2d); % Elevation
    
    % Fit Trigonometric States
    Cx(:,4) = ridgeLS(FBN, cos(psi_g), reg);
    Cx(:,5) = ridgeLS(FBN, sin(psi_g), reg);
    Cx(:,6) = ridgeLS(FBN, cos(theta_g), reg);
    Cx(:,7) = ridgeLS(FBN, sin(theta_g), reg);
end

% D. Fit Controls
Cu(:,1) = ridgeLS(FBN, tan(phi), reg); 
if Nu > 1
    Cu(:,2) = 0;
    Cu(:,3) = 0;
end

% -------------------------------------------------------------------------
% 5. Output Plot
% -------------------------------------------------------------------------
if isfield(data, 'plotinit') && data.plotinit
    PlotInitialGuess(Cx, Cu, data, uavIdx);
end

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos, vel, psi, phi] = evaluateBezierPhysics(p0, psi0, p3, psi3, R, t, V, g, phi_lim)
    % Inputs p0, p3 must be 1x2 row vectors
    
    % Bezier Control Points
    p1 = p0 + R * [cos(psi0), sin(psi0)];
    p2 = p3 - R * [cos(psi3), sin(psi3)];
    
    nt = length(t);
    pos = zeros(nt, 2);
    d_pos = zeros(nt, 2);  
    dd_pos = zeros(nt, 2); 
    
    for i = 1:nt
        tau = t(i);
        mt  = 1 - tau;
        
        pos(i,:)    = (mt^3)*p0 + 3*mt^2*tau*p1 + 3*mt*tau^2*p2 + tau^3*p3;
        d_pos(i,:)  = 3*mt^2*(p1-p0) + 6*mt*tau*(p2-p1) + 3*tau^2*(p3-p2);
        dd_pos(i,:) = 6*mt*(p2-2*p1+p0) + 6*tau*(p3-2*p2+p1);
    end
    
    psi = atan2(d_pos(:,2), d_pos(:,1));
    
    num = d_pos(:,1).*dd_pos(:,2) - d_pos(:,2).*dd_pos(:,1);
    den = (d_pos(:,1).^2 + d_pos(:,2).^2).^(1.5);
    K   = num ./ max(den, 1e-6);
    
    phi = atan( (V^2 / g) * K );
    phi = max(min(phi, phi_lim), -phi_lim);
    
    vel = [V * cos(psi), V * sin(psi)];
end

function L = computeBezierLength(p0, psi0, p3, psi3, R)
    % Gauss-Legendre Quadrature (10 points)
    w = [0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996, 0.295524224714753, ...
         0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688];
    x = [-0.973906528517172, -0.865063366688985, -0.679409568299024, -0.433395394129247, -0.148874338981631, ...
          0.148874338981631,  0.433395394129247,  0.679409568299024,  0.865063366688985,  0.973906528517172];
      
    p1 = p0 + R * [cos(psi0), sin(psi0)];
    p2 = p3 - R * [cos(psi3), sin(psi3)];
    
    L = 0;
    for i = 1:10
        tau = 0.5 * x(i) + 0.5;
        mt  = 1 - tau;
        d_vec = 3*mt^2*(p1-p0) + 6*mt*tau*(p2-p1) + 3*tau^2*(p3-p2);
        L = L + w(i) * norm(d_vec);
    end
    L = L * 0.5;
end

function c = ridgeLS(A, b, reg)
    n = size(A,2);
    c = (A.'*A + reg*eye(n)) \ (A.'*b);
end

% --- PLOTTING FUNCTION ---

function PlotInitialGuess(Cx, Cu, data, uavIdx)
    % Reconstruct smooth path for plotting
    t_plot = linspace(0, 1, 200)';
    BN = bernsteinMatrix_local(data.N_B, t_plot);
    FN = FourierMatrix_local(data.N_F, t_plot);
    Basis = [BN, FN];
    
    X_recon = Basis * Cx;
    U_recon = Basis * Cu;
    
    X_A     = X_recon(:,1);
    Y_A     = X_recon(:,2);
    Psi_A   = X_recon(:,3);
    tanPhiA = U_recon(:,1);
    
    idx_str = num2str(uavIdx);
    uavName = ['UAV' idx_str];
    Z_A     = data.(['z_A' idx_str]) * ones(size(X_A));
    
    % Pack temp traj struct for InFOV
    traj = struct();
    traj.(uavName).X_A     = X_A; 
    traj.(uavName).Y_A     = Y_A;
    traj.(uavName).Z_A     = Z_A;
    traj.(uavName).Psi_A   = Psi_A;
    traj.(uavName).Phi_A   = atan(tanPhiA);
    traj.(uavName).tan_Phi_A = tanPhiA;
    
    if data.isGimballed
        traj.(uavName).cPsi_g_A   = X_recon(:,4);
        traj.(uavName).sPsi_g_A   = X_recon(:,5);
        traj.(uavName).cTheta_g_A = X_recon(:,6);
        traj.(uavName).sTheta_g_A = X_recon(:,7);
    else
        % Dummy gimbal values if InFOV expects them
        traj.(uavName).cPsi_g_A   = zeros(size(X_A));
        traj.(uavName).sPsi_g_A   = zeros(size(X_A));
        traj.(uavName).cTheta_g_A = zeros(size(X_A));
        traj.(uavName).sTheta_g_A = zeros(size(X_A));
    end

    % Get Targets (Interpolated to plot grid)
    targets = [data.x_T(1,end), data.x_T(2,end), 0;  
               data.X_B(1,end), data.X_B(2,end), 0];
    
    % Check FOV
    camPos = [X_A.'; Y_A.'; Z_A.'];
    visMat = zeros(2, size(X_A,1));
    visMat(1,:) = InFOV(camPos, targets(1,:)', traj, uavIdx, data) > 1e-2;
    visMat(2,:) = InFOV(camPos, targets(2,:)', traj, uavIdx, data) > 1e-2;
    
    % Coloring
    baseColors = lines(2);
    posColors = zeros(size(X_A,1), 3);
    for k = 1:size(X_A,1)
        whichTargets = find(visMat(:,k));
        if isempty(whichTargets)
            posColors(k,:) = [0 0 0]; % Black = Empty
        else
            posColors(k,:) = mean(baseColors(whichTargets,:),1);
        end
    end
    
    % Actual Plot
    figName = sprintf('Initial Guess (Feasible) - UAV %d', uavIdx);
    figure('Name', figName, 'NumberTitle','off');
    hold on; grid on; axis equal;
    
    % Plot Trajectory Points
    scatter(X_A, Y_A, 20, posColors, 'filled');
    
    % Plot Targets
    plot(targets(1,1), targets(1,2), 'p', 'MarkerSize',12, 'MarkerFaceColor',baseColors(1,:), 'DisplayName','Target');
    plot(targets(2,1), targets(2,2), 's', 'MarkerSize',10, 'MarkerFaceColor',baseColors(2,:), 'DisplayName','Boat');
    
    % Plot Constraints
    theta = linspace(0,2*pi,100);
    r_NFZ = data.r_NFZ;
    rCom  = data.rCom;
    
    % NFZ around Target
    plot(targets(1,1)+r_NFZ*cos(theta), targets(1,2)+r_NFZ*sin(theta), '--k', 'LineWidth',1.5, 'DisplayName','NFZ');
    
    % Comm Range around Boat
    plot(targets(2,1)+rCom*cos(theta), targets(2,2)+rCom*sin(theta), '-.g', 'LineWidth',1.5, 'DisplayName','Comm Range');
    
    legend('Location','bestoutside');
    xlabel('X [m]'); ylabel('Y [m]');
    title(figName);
    hold off;
end

% --- BASIS DEFINITIONS (Must match Main Script) ---

function BN = bernsteinMatrix_local(N, t)
    t = t(:);
    BN = zeros(length(t), N + 1);
    for k = 0:N
        log_coeff = gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1);
        val = exp(log_coeff + k.*log(t + 1e-16) + (N-k).*log(1 - t + 1e-16));
        BN(:, k + 1) = val;
    end
    BN(1,1)=1; BN(1,2:end)=0;
    BN(end,end)=1; BN(end,1:end-1)=0;
end

function FN = FourierMatrix_local(N, t)
    t = t(:);
    w_base = 2*pi*t;
    FN = zeros(length(t), 2*N);
    for k = 1:N
        w = w_base * k;
        FN(:, (k-1)*2+1) = cos(w);
        FN(:, (k-1)*2+2) = sin(w);
    end
end