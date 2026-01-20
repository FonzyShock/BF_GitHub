function InFOV_ind = InFOV(cam_pos, cam_target, traj, uavIdx, data)
% Pure MATLAB version (No CasADi)

[~, N] = size(cam_pos);
if size(cam_target,2)==1
    cam_target  = repmat(cam_target,1,N);
end

idx_str = num2str(uavIdx);
     
uavStr  = ['UAV' idx_str];
psi_col = traj.(uavStr).Psi_A(:);    
phi_col = traj.(uavStr).Phi_A(:);    

if size(psi_col,1) ~= size(cam_pos,2)
    psi_col = data.BN * psi_col; 
    phi_col = data.BN * phi_col;
end

psi  = psi_col.';   
phi  = phi_col.';
cPsi = cos(psi);  
sPsi = sin(psi);
cPhi = cos(phi);  
sPhi = sin(phi);

% Relative in NED frame
rel_ned = cam_target - cam_pos;
x = rel_ned(1,:);  
y = rel_ned(2,:);  
z = rel_ned(3,:);

% Rotate into body frame
x1 = cPsi.*x + sPsi.*y;
y1 = -sPsi.*x + cPsi.*y;
z1 = z;

Xb = x1;
Yb = cPhi.*y1 + sPhi.*z1;
Zb = -sPhi.*y1 + cPhi.*z1;

rel_b = vertcat(Xb, Yb, Zb); 

if data.isGimballed
    cPsi_g = traj.(uavStr).cPsi_g_A.';
    sPsi_g = traj.(uavStr).sPsi_g_A.';  
    cTheta_g = traj.(uavStr).cTheta_g_A.'; 
    sTheta_g = traj.(uavStr).sTheta_g_A.'; 

    if size(cPsi_g,2) ~= size(cam_pos,2)
        cPsi_g   = data.BN * cPsi_g;  
        sPsi_g   = data.BN * sPsi_g;
        cTheta_g = data.BN * cTheta_g;  
        sTheta_g = data.BN * sTheta_g;
    end

    zero_elem = zeros(size(cPsi_g));

    R11 =  cPsi_g.*cTheta_g; R12 = -sPsi_g;        R13 =  cPsi_g.*sTheta_g;
    R21 =  sPsi_g.*cTheta_g; R22 =  cPsi_g;        R23 =  sPsi_g.*sTheta_g;
    R31 =  -sTheta_g;        R32 =  zero_elem;     R33 =  cTheta_g;

    xg = R11.*Xb + R12.*Yb + R13.*Zb;
    yg = R21.*Xb + R22.*Yb + R23.*Zb;
    zg = R31.*Xb + R32.*Yb + R33.*Zb;

    rel_b = vertcat(xg, yg, zg);
end

rel_az = atan2(rel_b(2,:), rel_b(1,:));
rel_el = atan2(-rel_b(3,:), sqrt(rel_b(1,:).^2 + rel_b(2,:).^2));

% Sigmoid FOV approximation
k         = 10;
EXP_MIN   = -700;  
EXP_MAX   = 700;
lim_az    = data.az_lim;
lim_el_up = data.el_lim_up;
lim_el_dn = data.el_lim_dn;

clip   = @(v) min(max(v, EXP_MIN), EXP_MAX);
sig    = @(x, thr) 1./(1 + exp(clip( k*(x - thr))));

el_up  = sig(rel_el, lim_el_up);
el_dn  = 1 - sig(rel_el, lim_el_dn);
el_FOV = el_up .* el_dn;

az_FOV = sig(abs(rel_az), lim_az/2);

InFOV_ind = (el_FOV .* az_FOV).';  
end