function [el_lim, az_lim, ang_res] = calc_FOV(boresight_elev, res_x, res_y, focal_length_mm, pixel_pitch_um)
    % calc_FOV calculates horizontal/vertical FOV, elevation limits, and angular resolution
    % for a given camera configuration in NED frame.
    %
    % INPUTS:
    %   boresight_elev    - boresight elevation angle in degrees (NED frame)
    %   res_x             - image width in pixels (default: 2592)
    %   res_y             - image height in pixels (default: 1944)
    %   focal_length_mm   - focal length in mm (default: 6.0 mm)
    %   pixel_pitch_um    - pixel pitch in microns (default: 2.2 µm)
    %
    % OUTPUTS:
    %   el_lim       - struct with elevation limits (deg, NED frame):
    %                    .max (most upward), .min (most downward)
    %   az_lim          - horizontal field of view in degrees
    %   ang_res           - struct with angular resolution (deg/pixel)

    % Set defaults if missing
   if nargin < 1 || isempty(boresight_elev)
        boresight_elev = 0;
    end
    if nargin < 2 || isempty(res_x)
        res_x = 2592;
    end
    if nargin < 3 || isempty(res_y)
        res_y = 1944;
    end
    if nargin < 4 || isempty(focal_length_mm)
        focal_length_mm = 6.0;
    end
    if nargin < 5 || isempty(pixel_pitch_um)
        pixel_pitch_um = 2.2;
    end

    % Convert pixel pitch to mm
    pixel_pitch_mm = pixel_pitch_um / 1000;

    % Sensor physical dimensions in mm
    sensor_width_mm = res_x * pixel_pitch_mm;
    sensor_height_mm = res_y * pixel_pitch_mm;

    % Compute horizontal and vertical FOV (in degrees)
    theta_h = 2 * atand(sensor_width_mm / (2 * focal_length_mm));
    theta_v = 2 * atand(sensor_height_mm / (2 * focal_length_mm));

    % Round FOVs to 1 decimal place
    theta_h = round(theta_h, 1);
    theta_v = round(theta_v, 1);

    % Angular resolution (degrees per pixel)
    ang_res.horizontal = round(theta_h / res_x, 4);
    ang_res.vertical   = round(theta_v / res_y, 4);

    % Elevation angle limits in degrees (relative to NED horizontal)
    el_lim.up = round(boresight_elev + theta_v / 2, 1);
    el_lim.dn = round(boresight_elev - theta_v / 2, 1);
    az_lim = theta_h; % Don't cut in half for inFOV.m code use. It does it there!
end
