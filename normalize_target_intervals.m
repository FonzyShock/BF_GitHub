function out = normalize_target_intervals(data)
% Normalize target intervals for 3D trajectories (boat + targets)

num_targets = data.num_targets; % includes boat already
t = data.num_time_steps;
time_opt = data.time_opt;
tf = data.tf;
sample_rate = data.V_A / data.dps; % Hz
sample_freq = 1 / sample_rate;     % seconds

% Extract Boat (USV) and Target trajectories
x_B = data.X_B(1,:);
y_B = data.X_B(2,:);
z_B = data.X_B(3,:);

x_T_all = squeeze(data.x_T(1,:,:)); % [N x num_targets-1] ideally
y_T_all = squeeze(data.x_T(2,:,:));
z_T_all = squeeze(data.x_T(3,:,:));

if isvector(x_T_all)
    % Handle 1-target case
    x_T_all = x_T_all(:);
    y_T_all = y_T_all(:);
    z_T_all = z_T_all(:);
end

time_steps = t;
uav_time = 0:sample_freq:(time_steps-1)*sample_freq;
target_time_range = time_opt * tf;

% Interpolate Boat
x_B_interp = interp1(target_time_range, x_B, uav_time, 'linear', 'extrap');
y_B_interp = interp1(target_time_range, y_B, uav_time, 'linear', 'extrap');
z_B_interp = interp1(target_time_range, z_B, uav_time, 'linear', 'extrap');

% Interpolate Targets
x_T_interp = zeros(num_targets-1, length(uav_time));
y_T_interp = zeros(num_targets-1, length(uav_time));
z_T_interp = zeros(num_targets-1, length(uav_time));

for tgt_idx = 1:num_targets-1
    x_T_interp(tgt_idx,:) = interp1(target_time_range, x_T_all(:,tgt_idx), uav_time, 'linear', 'extrap');
    y_T_interp(tgt_idx,:) = interp1(target_time_range, y_T_all(:,tgt_idx), uav_time, 'linear', 'extrap');
    z_T_interp(tgt_idx,:) = interp1(target_time_range, z_T_all(:,tgt_idx), uav_time, 'linear', 'extrap');
end

% Check and clear NaNs
arrays = {x_B_interp, y_B_interp, z_B_interp, x_T_interp, y_T_interp, z_T_interp};

for i = 1:numel(arrays)
    current_array = arrays{i};
    for tgt_or_row = 1:size(current_array,1)
        for j = 1:size(current_array,2)
            if isnan(current_array(tgt_or_row,j))
                if j == 1
                    current_array(tgt_or_row,j) = current_array(tgt_or_row,find(~isnan(current_array(tgt_or_row,:)), 1, 'first'));
                else
                    current_array(tgt_or_row,j) = current_array(tgt_or_row,j-1);
                end
            end
        end
    end
    arrays{i} = current_array;
end

% Unpack arrays
x_B_interp = arrays{1};
y_B_interp = arrays{2};
z_B_interp = arrays{3};
x_T_interp = arrays{4};
y_T_interp = arrays{5};
z_T_interp = arrays{6};

% Preallocate output array
out = zeros(num_targets, 3, length(uav_time));

% Boat 
for t_idx = 1:length(uav_time)
    out(1,:,t_idx) = [x_B_interp(t_idx), y_B_interp(t_idx), z_B_interp(t_idx)];
end

% targets
for tgt_idx = 1:(num_targets-1)
    for t_idx = 1:length(uav_time)
        out(tgt_idx+1,:,t_idx) = [x_T_interp(tgt_idx,t_idx), ...
                                  y_T_interp(tgt_idx,t_idx), ...
                                  z_T_interp(tgt_idx,t_idx)];
    end
end

end
