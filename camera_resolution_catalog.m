function CAM = camera_resolution_catalog()
% CAMERA_RESOLUTION_CATALOG
% Returns an array of structs with fields:
%   .name, .res_x, .res_y, .total_MP

    CAM = struct('name',{},'res_x',{},'res_y',{},'total_MP',{});

% % -------- Legacy / Baseline (4:3 & 5:4) --------
% CAM(end+1) = struct('name','VGA_640x480',   'res_x', 640, 'res_y', 480,  'total_MP', (640*480)/1e6);    % ~0.31 MP
% CAM(end+1) = struct('name','SVGA_800x600',  'res_x', 800, 'res_y', 600,  'total_MP', (800*600)/1e6);    % ~0.48 MP
% CAM(end+1) = struct('name','XGA_1024x768',  'res_x',1024, 'res_y', 768,  'total_MP', (1024*768)/1e6);   % ~0.79 MP
% CAM(end+1) = struct('name','SXGA_1280x1024','res_x',1280,'res_y',1024, 'total_MP', (1280*1024)/1e6);   % ~1.31 MP
% CAM(end+1) = struct('name','UXGA_1600x1200','res_x',1600,'res_y',1200, 'total_MP', (1600*1200)/1e6);   % ~1.92 MP
% CAM(end+1) = struct('name','QXGA_2048x1536','res_x',2048,'res_y',1536, 'total_MP', (2048*1536)/1e6);   % ~3.15 MP

% -------- 16:9 Video Family --------
CAM(end+1) = struct('name','HDplus_1600x900',   'res_x',1600,'res_y', 900, 'total_MP', (1600*900)/1e6);    % ~1.44 MP
CAM(end+1) = struct('name','QHD_2560x1440',     'res_x',2560,'res_y',1440,'total_MP', (2560*1440)/1e6);    % ~3.69 MP
CAM(end+1) = struct('name','WQXGA_2560x1600',   'res_x',2560,'res_y',1600,'total_MP', (2560*1600)/1e6);    % ~4.10 MP (16:10)
CAM(end+1) = struct('name','4K_UHD_3840x2160',  'res_x',3840,'res_y',2160,'total_MP', (3840*2160)/1e6);    % ~8.29 MP
CAM(end+1) = struct('name','DCI_4K_4096x2160',  'res_x',4096,'res_y',2160,'total_MP', (4096*2160)/1e6);    % ~8.85 MP
CAM(end+1) = struct('name','5K_5120x2880',      'res_x',5120,'res_y',2880,'total_MP', (5120*2880)/1e6);    % ~14.75 MP
CAM(end+1) = struct('name','8K_UHD_7680x4320',  'res_x',7680,'res_y',4320,'total_MP', (7680*4320)/1e6);    % ~33.18 MP
% CAM(end+1) = struct('name','8K_DCI_8192x4320',  'res_x',8192,'res_y',4320,'total_MP', (8192*4320)/1e6);    % ~35.39 MP

% % -------- Common Still-Photo 4:3 Sensors --------
% CAM(end+1) = struct('name','5MP_2592x1944',   'res_x',2592,'res_y',1944,'total_MP', (2592*1944)/1e6);    % ~5.04 MP
% CAM(end+1) = struct('name','8MP_3264x2448',   'res_x',3264,'res_y',2448,'total_MP', (3264*2448)/1e6);    % ~7.99 MP
% CAM(end+1) = struct('name','12MP_4000x3000',  'res_x',4000,'res_y',3000,'total_MP', (4000*3000)/1e6);    % 12.0 MP
% CAM(end+1) = struct('name','16MP_4608x3456',  'res_x',4608,'res_y',3456,'total_MP', (4608*3456)/1e6);    % ~15.93 MP
% CAM(end+1) = struct('name','24MP_6000x4000',  'res_x',6000,'res_y',4000,'total_MP', (6000*4000)/1e6);    % 24.0 MP
% CAM(end+1) = struct('name','32MP_6560x4928',  'res_x',6560,'res_y',4928,'total_MP', (6560*4928)/1e6);    % ~32.34 MP


    % Optional warning if MP and pixels diverge materially
    for k = 1:numel(CAM)
        approx = (double(CAM(k).res_x)*double(CAM(k).res_y))/1e6;
        if abs(approx - CAM(k).total_MP) > 0.5
            warning('CAM %s: total_MP=%.2f vs computed=%.2f MP', CAM(k).name, CAM(k).total_MP, approx);
        end
    end
end
