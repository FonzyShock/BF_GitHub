function [tr_crlb_AB, tr_crlb_AT, data] = calc_CRLB_coop(data, results)
% CALC_CRLB_COOP_CASADI  Numeric CRLB print‑out for cooperative multi‑UAV
%   [tr_crlb_AB, tr_crlb_AT, data] = calc_CRLB_coop_casadi(data, results)
%
% Inputs:
%   data    – problem struct (must include .num_uav, .tf, etc.)
%   results – struct with field xout (decision vector) and jout (cost)
%
% Outputs:
%   tr_crlb_AB – combined trace CRLB to Boat
%   tr_crlb_AT – combined trace CRLB to Target
%   data       – augmented with per‑UAV and combined CRLB fields

import casadi.*

% Initialize combined FIM accumulators
num_uav = data.num_uav;
FIM_AT_sum = zeros(3,3);
FIM_AB_sum = zeros(3,3);
tol = data.tol;

try
    % Reconstruct trajectories (numeric) via CasADi helper
    X   = results.xout;
    traj = extractTrajectories_Cooperative(X, data);

    fprintf('Cooperative UAV Data\n');
    % Per‑UAV CRLB
    for uavIdx = 1:data.num_uav
        idx = num2str(uavIdx);
        uavField = ['UAV' idx];
        % Skip missing
        if ~isfield(traj, uavField)
            fprintf('  [UAV %s] Trajectory missing; skipping CRLB.\n', idx);
            continue;
        end

        % --- Compute FIM ---
        [FIM_AT, FIM_AB, tr_AT, tr_AB] = calc_FIM_CRLB(traj, idx, data);
        FIM_AT = full(FIM_AT);
        FIM_AB = full(FIM_AB);

        % % --- CRLB via safe pseudo-inverse ---
        % crlb_AT = safeInverse(FIM_AT, tol);
        % crlb_AB = safeInverse(FIM_AB, tol);

        % % Trace of CRLB
        % tr_AT = trace(crlb_AT);
        % tr_AB = trace(crlb_AB);

         % Accumulate for combined result
    FIM_AT_sum = FIM_AT_sum + full(FIM_AT);
    FIM_AB_sum = FIM_AB_sum + full(FIM_AB);

        % Store in data
        data.(['tr_crlb_AT_coop_UAV' idx]) = tr_AT;
        data.(['tr_crlb_AB_coop_UAV' idx]) = tr_AB;

         % Print per‑UAV CRLBs
    fprintf('  UAV %s - Target CRLB trace: %2.6f\n', idx, double(full(evalf(tr_AT))));
    fprintf('  UAV %s - Boat   CRLB trace: %2.6f\n', idx, double(full(evalf(tr_AB))));
    end

    % Combined CRLB traces
    % Combined CRLBs
crlb_AT_comb = safeInverse(FIM_AT_sum, tol);
crlb_AB_comb = safeInverse(FIM_AB_sum, tol);

 tr_crlb_AT = trace(crlb_AT_comb);
    tr_crlb_AB = trace(crlb_AB_comb);

    tr_crlb_AT = trace(inv(FIM_AT_sum));
    tr_crlb_AB = trace(inv(FIM_AB_sum));

    % Store combined in data
    data.tr_crlb_AT_coop_combined = tr_crlb_AT;
    data.tr_crlb_AB_coop_combined = tr_crlb_AB;

    % Print combined sqrt‐PCRLB traces and summary
    fprintf('  sqrt(PCRLB trace) to Target ‑ Combined: %2.6f\n', sqrt(double(full(evalf(tr_crlb_AT)))));
    fprintf('  sqrt(PCRLB trace) to Boat   ‑ Combined: %2.6f\n', sqrt(double(full(evalf(tr_crlb_AB)))));
    fprintf('  Optimization Cost: %2.4f\n', results.jout);
    total_time = data.tf * (1 + X(end));
    fprintf('  Mission Time: %2.2f s\n', total_time);

catch ME
    fprintf('Error computing cooperative CRLBs: %s\n', ME.message);
end
end
