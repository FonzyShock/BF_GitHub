function traj = extractTrajectories_Cooperative(X, data)
% Pure MATLAB version (No CasADi)

num_uav     = data.num_uav;
N_B         = data.N_B;
N_tot       = data.N_tot;
Nx          = data.Nx;
Nu          = data.Nu;
isGimballed = data.isGimballed;

nCx_per_uav = Nx * N_tot;
nCu_per_uav = Nu * N_tot;

start_Cx = 1;
end_Cx   = num_uav * nCx_per_uav;
start_Cu = end_Cx + 1;
end_Cu   = end_Cx + num_uav * nCu_per_uav;

Cx_all = X(start_Cx : end_Cx);
Cu_all = X(start_Cu : end_Cu);

Cx_mat = reshape(Cx_all, N_tot, num_uav * Nx);
Cu_mat = reshape(Cu_all, N_tot, num_uav * Nu);

BN  = data.BN;
FN  = data.FN;
DBm = data.DBm;
DFm = data.DFm;

traj = struct();

for i = 1:num_uav
    uav_name = sprintf('UAV%d', i);
    
    s_idx = (i-1)*Nx + (1:Nx);
    c_idx = (i-1)*Nu + (1:Nu);
    
    Cxi = Cx_mat(:, s_idx);
    Cui = Cu_mat(:, c_idx);
    
    CxB = Cxi(1 : N_B+1, :);
    CxF = Cxi(N_B+2 : end, :);
    
    CuB = Cui(1 : N_B+1, :);
    CuF = Cui(N_B+2 : end, :);
    
    X_recon = BN * CxB + FN * CxF;
    U_recon = BN * CuB + FN * CuF;
    dX_recon = DBm * CxB + DFm * CxF;
    
    traj.(uav_name).X_A     = X_recon(:,1);
    traj.(uav_name).Y_A     = X_recon(:,2);
    traj.(uav_name).Psi_A   = X_recon(:,3);
    
    z_field = sprintf('z_A%d', i);
    traj.(uav_name).Z_A     = data.(z_field) * ones(size(BN,1),1);
    
    traj.(uav_name).dX_A    = dX_recon(:,1);
    traj.(uav_name).dY_A    = dX_recon(:,2);
    traj.(uav_name).dPsi_A  = dX_recon(:,3);
    
    traj.(uav_name).tan_Phi_A = U_recon(:,1);
    traj.(uav_name).Phi_A     = atan(U_recon(:,1));
    
    if isGimballed
        traj.(uav_name).cPsi_g_A    = X_recon(:,4);
        traj.(uav_name).sPsi_g_A    = X_recon(:,5);
        traj.(uav_name).cTheta_g_A  = X_recon(:,6);
        traj.(uav_name).sTheta_g_A  = X_recon(:,7);

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
end