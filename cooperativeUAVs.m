function [xOut, Jout, data_out] = cooperativeUAVs(data, results)
% FMINCON Version

data_out = data;
num_uav  = data.num_uav;
N_tot    = data.N_tot;
Nx       = data.Nx;
Nu       = data.Nu;
MissionAreaBound = data.MissionAreaBound;

% Build Initial Guess
Cx_all = []; Cu_all = []; alpha = 0; X_A_0 = []; u_0 = [];
size_Cx_i = Nx * N_tot;
size_Cu_i = Nu * N_tot;

for uavIdx = 1:num_uav
    idx = num2str(uavIdx);
    xi_name = ['xOut' idx];
    
    if isfield(results, xi_name) && ~isempty(results.(xi_name))
        Xi = results.(xi_name);
        Ci = Xi(1 : size_Cx_i);
        Ui = Xi(size_Cx_i+1 : size_Cx_i+size_Cu_i);
        alpha = alpha + Xi(end);
    else
        [Ci0, Ui0] = initPath(data, uavIdx);
        Ci = reshape(Ci0, [], 1);
        Ui = reshape(Ui0, [], 1);
    end
    
    Cx_all = [Cx_all; Ci];
    Cu_all = [Cu_all; Ui];
    x0_A = data.(['x0_A' idx]);
    X_A_0 = [X_A_0; x0_A(:)];
    u_0 = [u_0; Ui(1)]; 
end

alpha_avg = alpha / num_uav;
X0 = [Cx_all; Cu_all; alpha_avg];
nVar = numel(X0);

% Bounds
x_lb = -inf(num_uav * Nx * N_tot, 1);
x_ub =  inf(num_uav * Nx * N_tot, 1);

for i = 1:num_uav
    offset = (i-1) * Nx * N_tot;
    x_lb(offset + 1 : offset + N_tot) = MissionAreaBound(1);
    x_ub(offset + 1 : offset + N_tot) = MissionAreaBound(2);
    x_lb(offset + N_tot + 1 : offset + 2*N_tot) = MissionAreaBound(3);
    x_ub(offset + N_tot + 1 : offset + 2*N_tot) = MissionAreaBound(4);
    x_lb(offset + 2*N_tot + 1 : offset + 3*N_tot) = -pi;
    x_ub(offset + 2*N_tot + 1 : offset + 3*N_tot) =  pi;
end

u_lb = -ones(num_uav * Nu * N_tot, 1) * data.control_limit;
u_ub = -u_lb;
lb = [x_lb; u_lb; data.alpha_min];
ub = [x_ub; u_ub; data.alpha_max];

% Equality Constraints
nEq = num_uav * (Nx + Nu);
Aeq = zeros(nEq, nVar);
beq = zeros(nEq, 1);

Basis_t0 = [data.BN(1,:), data.FN(1,:)]; 
row_count = 0;

for i = 1:num_uav
    for s = 1:Nx
        row_count = row_count + 1;
        block_idx = (i-1)*Nx + (s-1);
        col_start = block_idx * N_tot + 1;
        col_end   = (block_idx + 1) * N_tot;
        Aeq(row_count, col_start:col_end) = Basis_t0;
        
        x0 = data.(['x0_A' num2str(i)]);
        if s <= 3, beq(row_count) = x0(s);
        else, if s==4 || s==6, beq(row_count)=1; else, beq(row_count)=0; end
        end
    end
end

state_offset = num_uav * Nx * N_tot;
for i = 1:num_uav
    for c = 1:Nu
        row_count = row_count + 1;
        block_idx = (i-1)*Nu + (c-1);
        col_start = state_offset + block_idx * N_tot + 1;
        col_end   = state_offset + (block_idx + 1) * N_tot;
        Aeq(row_count, col_start:col_end) = Basis_t0;
        beq(row_count) = X0(col_start); 
    end
end

% Options
options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ... 
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 2e5, ...
    'MaxIterations', 3000, ...
    'ConstraintTolerance', 1e-4);

fun = @(x) costFunc(x, data);
nonlcon = @(x) nonlcon_coop_uav(x, data);

try
    [xOut, Jout, exitflag, output] = fmincon(fun, X0, [], [], Aeq, beq, lb, ub, nonlcon, options);
    data.convergeFlag = (exitflag > 0);
    data.output.status = output.message;
catch
    xOut = []; Jout = [];
    data.convergeFlag = 0;
    data.output.status = "Failed";
end
data_out = data;
end

function J = costFunc(X, data)
traj = extractTrajectories_Cooperative(X, data);
SUM_AT = zeros(3,3);
SUM_AB = zeros(3,3);

for i = 1:data.num_uav
    [F_AT, F_AB] = calc_FIM_CRLB(traj, num2str(i), data);
    SUM_AT = SUM_AT + F_AT;
    SUM_AB = SUM_AB + F_AB;
end

eps_fac = 1e-6;
J_AT = trace(inv(SUM_AT + eps_fac*eye(3)));
J_AB = trace(inv(SUM_AB + eps_fac*eye(3)));

w = data.priority;
J = w * J_AT + (1-w)*J_AB + data.time_weight * X(end);
end