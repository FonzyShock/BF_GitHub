function traj = extractTrajectories_nonCooperative(X, data)
% Pure MATLAB version (No CasADi)

N_B         = data.N_B;
N_tot       = data.N_tot;        
Nx          = data.Nx;      
Nu          = data.Nu;      
isGimballed = data.isGimballed;
idx_num     = data.idx;     
uav_name    = ['UAV' num2str(idx_num)];  

% --- 1. Extract State Coefficients ---
nCx    = Nx * N_tot;
Cx_vec = X(1:nCx);
Cx     = reshape(Cx_vec, N_tot, Nx);

CxB = Cx(1 : N_B+1, :);
CxF = Cx(N_B+2 : end, :);

% --- 2. Extract Control Coefficients ---
if isGimballed
    Cu_start = nCx + 1;
    Cu_end   = Cu_start + Nu*N_tot - 1;
else
    Cu_start = nCx + 1;
    Cu_end   = numel(X) - 1; % Exclude alpha
end

Cu_vec = X(Cu_start : Cu_end);
Cu     = reshape(Cu_vec, N_tot, Nu); 

CuB = Cu(1 : N_B+1, :);
CuF = Cu(N_B+2 : end, :);

% --- 3. Reconstruct Trajectories ---
BN = data.BN;  
FN = data.FN;  

X_recon = BN * CxB + FN * CxF;

traj = struct();
traj.(uav_name).X_A   = X_recon(:,1);
traj.(uav_name).Y_A   = X_recon(:,2);
traj.(uav_name).Psi_A = X_recon(:,3);

z_field = ['z_A' num2str(idx_num)];            
traj.(uav_name).Z_A = data.(z_field) * ones(size(BN,1), 1);

% --- 4. Reconstruct Derivatives ---
DBm = data.DBm; 
DFm = data.DFm;

dX_recon = DBm * CxB + DFm * CxF;

traj.(uav_name).dX_A   = dX_recon(:,1);
traj.(uav_name).dY_A   = dX_recon(:,2);
traj.(uav_name).dPsi_A = dX_recon(:,3);

% --- 5. Reconstruct Controls ---
U_recon = BN * CuB + FN * CuF;

traj.(uav_name).tan_Phi_A = U_recon(:,1);
traj.(uav_name).Phi_A     = atan(U_recon(:,1));

% --- 6. Gimballed-UAV Additional States ---
if isGimballed
    traj.(uav_name).cPsi_g_A   = X_recon(:,4);
    traj.(uav_name).sPsi_g_A   = X_recon(:,5);
    traj.(uav_name).cTheta_g_A = X_recon(:,6);
    traj.(uav_name).sTheta_g_A = X_recon(:,7);

    traj.(uav_name).dcPsi_g_A   = dX_recon(:,4);
    traj.(uav_name).dsPsi_g_A   = dX_recon(:,5);
    traj.(uav_name).dcTheta_g_A = dX_recon(:,6);
    traj.(uav_name).dsTheta_g_A = dX_recon(:,7);

    traj.(uav_name).u_Psi_g_A   = U_recon(:,2);
    traj.(uav_name).u_Theta_g_A = U_recon(:,3);

    traj.(uav_name).Psi_g_A   = atan2(X_recon(:,5), X_recon(:,4));
    traj.(uav_name).Theta_g_A = atan2(X_recon(:,7), X_recon(:,6));
end
end