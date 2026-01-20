function [c, ceq] = nonlcon_coop_uav(X, data)
% Pure MATLAB version (No CasADi)

V_A        = data.V_A;
x_T        = data.x_T;
X_B        = data.X_B;
rCom       = data.rCom;
r_NFZ      = data.r_NFZ;
g          = data.g;
num_uav    = data.num_uav;
isGimballed = data.isGimballed;

alpha = X(end);
tf    = data.tf * (1 + alpha);

traj = extractTrajectories_Cooperative(X, data);

c_all   = [];
ceq_all = [];

for i = 1:num_uav
    uavField = sprintf('UAV%d', i);
    u = traj.(uavField);

    ceq_all = [ceq_all;
        u.dPsi_A - (tf * g/V_A) .* u.tan_Phi_A;
        u.dX_A   - tf * V_A .* cos(u.Psi_A);
        u.dY_A   - tf * V_A .* sin(u.Psi_A)
    ];

    dx_AT = x_T(1) - u.X_A;
    dy_AT = x_T(2) - u.Y_A;
    NFZ_lim = r_NFZ^2 - (dx_AT.^2 + dy_AT.^2);

    dx_AB_end = X_B(1,end) - u.X_A(end);
    dy_AB_end = X_B(2,end) - u.Y_A(end);
    comm_end = (dx_AB_end^2 + dy_AB_end^2) - rCom^2;

    c_all = [c_all;
             comm_end;
             NFZ_lim;
    ];

    if isGimballed
        cPsi = u.cPsi_g_A;   
        sPsi = u.sPsi_g_A;
        dcPsi= u.dcPsi_g_A;  
        dsPsi= u.dsPsi_g_A;
        cTheta = u.cTheta_g_A; 
        sTheta = u.sTheta_g_A;
        dcTheta= u.dcTheta_g_A;
        dsTheta= u.dsTheta_g_A;
        uPsi_g = u.u_Psi_g_A;  
        uTheta_g = u.u_Theta_g_A;

        ceq_all = [ceq_all;
            dcPsi + sPsi   .* (tf * uPsi_g);
            dsPsi - cPsi   .* (tf * uPsi_g);
            dcTheta + sTheta   .* (tf * uTheta_g);
            dsTheta - cTheta   .* (tf * uTheta_g);
        ];
    end
end

c   = c_all;
ceq = ceq_all;
end