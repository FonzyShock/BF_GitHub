function [xOut, Jout, data_out] = nonCooperativeUAVs(data, uavIdx)
% FMINCON Version

data_out = data; 
N_tot    = data.N_tot; 
Nx       = data.Nx;
Nu       = data.Nu;
isGimbled           = data.isGimballed;
MissionAreaBound    = data.MissionAreaBound;
alpha_min           = data.alpha_min;
alpha_max           = data.alpha_max;

idx                 = num2str(uavIdx);
x0_A                = data.(['x0_A' idx]);

% Init Path using BF
[Cx0, Cu0]          = initPath(data, uavIdx);
Cx0                 = reshape(Cx0, Nx*N_tot, 1);
Cu0                 = reshape(Cu0, Nu*N_tot, 1);
alpha0              = 0; 

X0                  = [Cx0; Cu0; alpha0]; 
nVar                = (Nx+Nu)*N_tot + 1;

% Bounds
x_lb = -inf(Nx*N_tot, 1);
x_ub =  inf(Nx*N_tot, 1);
u_lb = -inf(Nu*N_tot, 1);
u_ub =  inf(Nu*N_tot, 1);

x_lb(1 : N_tot)             = MissionAreaBound(1); 
x_ub(1 : N_tot)             = MissionAreaBound(2);
x_lb(N_tot+1 : 2*N_tot)     = MissionAreaBound(3); 
x_ub(N_tot+1 : 2*N_tot)     = MissionAreaBound(4);
x_lb(2*N_tot+1 : 3*N_tot)   = -pi;                 
x_ub(2*N_tot+1 : 3*N_tot)   =  pi;

if isGimbled
    x_lb(3*N_tot+1 : end) = -1;
    x_ub(3*N_tot+1 : end) =  1;
    
    u_lb(1 : N_tot)           = -data.control_limit; 
    u_ub(1 : N_tot)           =  data.control_limit;
    u_lb(N_tot+1 : 2*N_tot)   =  data.u_psi_g_min;
    u_ub(N_tot+1 : 2*N_tot)   =  data.u_psi_g_max;
    u_lb(2*N_tot+1 : 3*N_tot) =  data.u_theta_g_min;
    u_ub(2*N_tot+1 : 3*N_tot) =  data.u_theta_g_max;
else
    u_lb = -ones(Nu*N_tot, 1) * data.control_limit;
    u_ub = -u_lb;
end

lb = [x_lb; u_lb; alpha_min];
ub = [x_ub; u_ub; alpha_max];

% Linear Equality (Initial Conditions)
Aeq = zeros(Nx+Nu, nVar);
Basis_t0 = [data.BN(1,:), data.FN(1,:)]; 

for i = 1:(Nx+Nu)
    col_start = (i-1)*N_tot + 1;
    col_end   = i*N_tot;
    Aeq(i, col_start:col_end) = Basis_t0;
end

if isGimbled
    cPsi_g0=cos(0); sPsi_g0=sin(0); cTheta_g0=cos(0); sTheta_g0=sin(0);
    beq = [x0_A(1:3).'; cPsi_g0; sPsi_g0; cTheta_g0; sTheta_g0; Cu0(1); 0; 0];
else
    beq = [x0_A(1:3).'; Cu0(1)]; 
end

% FMINCON Setup
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ... 
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 1e3, ...
    'MaxIterations', 3000, ...
    'ConstraintTolerance', 1e-3, ...
    'StepTolerance', 1e-8);

data.idx = idx;
fun = @(x) costFunc(x, data);
nonlcon = @(x) nonlcon_noncoop_uav(x, data);

try
    tic;
    [xOut, Jout, exitflag, output] = fmincon(fun, X0, [], [], Aeq, beq, lb, ub, nonlcon, options);
    data.output.solve_time = toc;
    
    data.convergeFlag = (exitflag > 0);
    data.output.status = output.message;
    data.output.success = (exitflag > 0);
    data.output.iterations = output.iterations;

catch ME
    warning('Optimization failed: %s', ME.message);
    xOut = []; Jout = [];
    data.convergeFlag = 0;
    data.output.status = "Failed";
end
end

function J = costFunc(X, data)
w   = data.priority;       
tw  = data.time_weight;    
traj = extractTrajectories_nonCooperative(X, data);
[FIM_AT, FIM_AB] = calc_FIM_CRLB(traj, data.idx, data);

eps_fac = 1e-6;
FIM_AT_reg = FIM_AT + eps_fac * norm(FIM_AT) * eye(3);
FIM_AB_reg = FIM_AB + eps_fac * norm(FIM_AB) * eye(3);

J_AT = trace(inv(FIM_AT_reg));
J_AB = trace(inv(FIM_AB_reg));

J = w * J_AT + (1 - w) * J_AB + tw * X(end);
end