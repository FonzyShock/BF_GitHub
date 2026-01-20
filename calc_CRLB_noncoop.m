function [crlb_AB_cell, crlb_AT_cell, data] = calc_CRLB_noncoop(data, results)
% CALC_CRLB_NONCOOP  Numeric CRLB print‑out for non‑cooperative UAVs
%   [crlb_AB_cell, crlb_AT_cell, data] = calc_CRLB_NONCOOP(data, results)
%
% Inputs:
%   data    – your problem struct (must include .num_uav, .tf, etc.)
%   results – struct with fields xOut1, Jout1, xOut2, Jout2, …  
%
% Outputs:
%   crlb_AB_cell – 1×num_uav cell of 3×3 numeric CRLB to Boat
%   crlb_AT_cell – 1×num_uav cell of 3×3 numeric CRLB to Target
%   data         – data struct augmented with fields:
%                  .crlb_AT_noncoop_UAVi, .crlb_AB_noncoop_UAVi,
%                  .tr_crlb_AT_noncoop_combined, .tr_crlb_AB_noncoop_combined

% Preallocate accumulators
tol = data.tol;
num_uav = data.num_uav;
FIM_AT_sum = zeros(3,3);
FIM_AB_sum = zeros(3,3);

crlb_AT_cell = cell(1,num_uav);
crlb_AB_cell = cell(1,num_uav);

fprintf('Non‑cooperative UAV Data\n');

for uavIdx = 1:num_uav
    idx = num2str(uavIdx);
    xfield = ['xOut' idx];
    jfield = ['Jout' idx];

    % Skip if missing
    if ~isfield(results, xfield) || isempty(results.(xfield))
        fprintf('  [UAV %s] No solution – skipping CRLB.\n', idx);
        continue;
    end

    % Numeric solution & cost
    Xc = results.(xfield);     % ( (Nx+Nu)*(N+1)+1 )×1
    Jc = results.(jfield);

    % Reconstruct trajectory numerically
    data.idx = idx;  % for the helper calls
    traj = extractTrajectories_nonCooperative(Xc, data);

    % Numeric FIM computation
    [FIM_AT, FIM_AB] = calc_FIM_CRLB(traj, idx, data);

    % Invert to get CRLB matrices
     % --- Handle singular or ill-conditioned FIMs robustly ---
    tol = 1e-10;
    crlb_AT = safeInverse(FIM_AT, tol);
    crlb_AB = safeInverse(FIM_AB, tol);

    crlb_AT_cell{uavIdx} = crlb_AT;
    crlb_AB_cell{uavIdx} = crlb_AB;

    % Trace (PCRLB cost term)
    tr_crlb_AT = trace(crlb_AT);
    tr_crlb_AB = trace(crlb_AB);

    % Accumulate for combined result
    FIM_AT_sum = FIM_AT_sum + full(FIM_AT);
    FIM_AB_sum = FIM_AB_sum + full(FIM_AB);

    % Store in data
    data.(['crlb_AT_noncoop_UAV' idx]) = crlb_AT;
    data.(['crlb_AB_noncoop_UAV' idx]) = crlb_AB;

    % Print per‑UAV CRLBs
    fprintf('  UAV %s - Target CRLB trace: %2.6f\n', idx, double(full(evalf(tr_crlb_AT))));
    fprintf('  UAV %s - Boat   CRLB trace: %2.6f\n', idx, double(full(evalf(tr_crlb_AB))));
end

% Combined CRLBs
crlb_AT_comb = safeInverse(FIM_AT_sum, tol);
crlb_AB_comb = safeInverse(FIM_AB_sum, tol);

tr_AT_comb   = trace(crlb_AT_comb);
tr_AB_comb   = trace(crlb_AB_comb);

data.tr_crlb_AT_noncoop_combined = tr_AT_comb;
data.tr_crlb_AB_noncoop_combined = tr_AB_comb;

% Final printout
fprintf('Combined Target sqrt(CRLB trace): %2.6f\n', sqrt(double(full(evalf(trace(crlb_AT_comb))))));
fprintf('Combined Boat   sqrt(CRLB trace): %2.6f\n', sqrt(double(full(evalf(trace(crlb_AB_comb))))));

% Also echo the last mission cost & time
if exist('Jc','var') && ~isempty(Jc)
    alpha = Xc(end);
    total_time = data.tf * (1 + alpha);
    fprintf('Optimization cost: %2.4f\n', Jc);
    fprintf('Mission time:     %2.2f s\n', total_time);
end
end
