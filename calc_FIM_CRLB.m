function [FIM_AT, FIM_AB, tr_crlb_AT, tr_crlb_AB, FIM_AT_9xN, FIM_AB_9xN] = calc_FIM_CRLB(traj, uavIdx, data)
% Pure MATLAB version (No CasADi)

X_T     = data.x_T;
X_B     = data.X_B;
Sigma   = data.Sigma;
lambda  = data.lambda;

if isfield(data,'P0_AT'), P0_AT = data.P0_AT; else, P0_AT = []; end
if isfield(data,'P0_AB'), P0_AB = data.P0_AB; else, P0_AB = []; end
if isfield(data, 'Q'), Q_AT = data.Q; Q_AB = data.Q; else, Q_AT = []; Q_AB = []; end

idx = sprintf("UAV%s", uavIdx);
if ~isfield(traj, idx), idx = uavIdx; end

x_A = traj.(idx).X_A;    
y_A = traj.(idx).Y_A;    
z_A = traj.(idx).Z_A;    

if ~isequal(length(X_B'), length(x_A))
    % If grids don't match, assume linear interp was handled or sizes are correct
    % For safety, one could add interp1 logic here for numeric arrays
end

cam_pos    = vertcat(x_A.', y_A.', z_A.'); 
cam_target = X_T; 
cam_boat   = X_B; 

AT_ind = InFOV(cam_pos, cam_target, traj, uavIdx, data);  
AB_ind = InFOV(cam_pos, cam_boat,   traj, uavIdx, data);  

[FIM_AT, FIM_AT_9xN] = compute3DFIM(cam_pos, cam_target, AT_ind, Sigma);
[FIM_AB, FIM_AB_9xN] = compute3DFIM(cam_pos, cam_boat,   AB_ind, Sigma);

F_prior_AT = zeros(3,3);
if ~isempty(P0_AT) && lambda ~= 0, F_prior_AT = lambda * inv(P0_AT); end

F_prior_AB = zeros(3,3);
if ~isempty(P0_AB) && lambda ~= 0, F_prior_AB = lambda * inv(P0_AB); end

N_steps = size(x_A, 1);
F_proc_AT = zeros(3,3);
if ~isempty(Q_AT), F_proc_AT = N_steps * inv(Q_AT); end

F_proc_AB = zeros(3,3);
if ~isempty(Q_AB), F_proc_AB = N_steps * inv(Q_AB); end

FIM_AT = F_prior_AT + FIM_AT; 
FIM_AB = F_prior_AB + FIM_AB; 

if ~isempty(FIM_AT_9xN)
    F_full_prior_AT = reshape(F_prior_AT + F_proc_AT, 9, 1);
    FIM_AT_9xN = F_full_prior_AT + FIM_AT_9xN;
end
if ~isempty(FIM_AB_9xN)
    F_full_prior_AB = reshape(F_prior_AB + F_proc_AB, 9, 1);
    FIM_AB_9xN = F_full_prior_AB + FIM_AB_9xN;
end

eps_fac     = 1e-6;
eps_min     = 1e-9;
epsilon_AT  = max(eps_min, eps_fac * norm(FIM_AT(:)));
epsilon_AB  = max(eps_min, eps_fac * norm(FIM_AB(:)));

FIM_AT_reg  = FIM_AT + epsilon_AT * eye(3);
FIM_AB_reg  = FIM_AB + epsilon_AB * eye(3);

tr_crlb_AT  = sum(diag(inv(FIM_AT_reg)));
tr_crlb_AB  = sum(diag(inv(FIM_AB_reg)));
end

function [FIM, FIM_cum9xN] = compute3DFIM(cam_pos, target_pos, ind, Sigma)
N = size(cam_pos, 2);
eps_guard = 1e-6; 

r           = target_pos - cam_pos;
r2          = sum(r.^2, 1) + eps_guard;
rxy2        = sum(r(1:2,:).^2, 1) + eps_guard;
rxy         = sqrt(rxy2);

H_az        = [ -r(2,:) ./ rxy2;                 
                 r(1,:) ./ rxy2;                 
                zeros(1,N) ];                 

H_el        = [ r(1,:).*r(3,:) ./ (r2 .* rxy);   
                r(2,:).*r(3,:) ./ (r2 .* rxy);   
                -rxy ./ r2 ];                    

w_az        = 1 / Sigma(1,1);
w_el        = 1 / Sigma(2,2);

s           = (sqrt(ind(:))).';

Ha          = H_az .* (ones(3,1) * s);
He          = H_el .* (ones(3,1) * s);

FIM         = w_az * (Ha * Ha.') + w_el * (He * He.');

FIM_cum9xN  = zeros(9, N);
F_cum       = zeros(3,3);

for k = 1:N
    Hak     = Ha(:,k);
    Hek     = He(:,k);
    Fk      = w_az * (Hak * Hak.') + w_el * (Hek * Hek.');  
    F_cum   = F_cum + Fk;                                   
    FIM_cum9xN(:,k) = reshape(F_cum, 9, 1);                 
end
end