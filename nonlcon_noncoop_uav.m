function [c, ceq] = nonlcon_noncoop_uav(X, data)
% Pure MATLAB version (No CasADi)

g        = data.g;
V_A      = data.V_A;
X_T      = data.x_T;
X_B      = data.X_B;
rCom     = data.rCom;
r_NFZ    = data.r_NFZ;
isGimballed = data.isGimballed;
alpha    = X(end);
tf       = data.tf * (1 + alpha);
area     = data.MissionAreaBound; 

idx = data.idx;             

traj = extractTrajectories_nonCooperative(X, data);
uav  = traj.(sprintf("UAV%s", idx));

x_A     = uav.X_A;      
y_A     = uav.Y_A;
psi_A   = uav.Psi_A;    
tanPhi  = uav.tan_Phi_A;
dX_A    = uav.dX_A;     
dY_A    = uav.dY_A;
dPsi_A  = uav.dPsi_A;

% Dynamics Equality Constraints
ceq = [ ...
    dPsi_A - tf*(g/V_A).*tanPhi;      
    dX_A   - tf*V_A.*cos(psi_A);      
    dY_A   - tf*V_A.*sin(psi_A)       
    ];

% Range Constraints
dx_AB_end = X_B(1,end) - x_A(end);
dy_AB_end = X_B(2,end) - y_A(end);
dH2_AB    = dx_AB_end.^2 + dy_AB_end.^2;

dx_AT = X_T(1,:) - x_A.';
dy_AT = X_T(2,:) - y_A.';
dH2_AT = dx_AT.^2 + dy_AT.^2;

c = [ ...
    dH2_AB - rCom^2;            
    (r_NFZ^2 - dH2_AT).';         
    ];

c = [c;
     area(1) - x_A;  
     x_A - area(2);  
     area(3) - y_A;
     y_A - area(4);
     tanPhi - data.control_limit; 
     -data.control_limit - tanPhi
    ];

if isGimballed
    cPsi        = uav.cPsi_g_A;  
    sPsi        = uav.sPsi_g_A;
    dcPsi       = uav.dcPsi_g_A; 
    dsPsi       = uav.dsPsi_g_A;
    cTheta      = uav.cTheta_g_A; 
    sTheta      = uav.sTheta_g_A;
    dcTheta     = uav.dcTheta_g_A; 
    dsTheta     = uav.dsTheta_g_A;
    uPsi_g      = uav.u_Psi_g_A;    
    uTheta_g    = uav.u_Theta_g_A;

    ceq = [ ceq;
        dcPsi + sPsi.*(tf*uPsi_g);    
        dsPsi - cPsi.*(tf*uPsi_g);    
        dcTheta + sTheta.*(tf*uTheta_g);
        dsTheta - cTheta.*(tf*uTheta_g);
        ];
end
end