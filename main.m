function main(varargin)
% main function
% Before run this function, 'preprocess.m' should run first 
% and four .mat file (nodes, mask, rot_axis, filter_s) should 
% exist in current current directory.
clear global

%% set config
% whether write cloudpoints data into txt
% format(each row): [x_idx y_idx z_idx]
config.write_into_txt = false;
% realtime display 3d cloudpoints result
config.realtime_disp = false;
% read all images into memory previous
config.read_img_prev = true;
% dynamic adjust size of ROI mask (can speed up)
config.mask_dynamc_adj = true;
% timing
config.timing = true;
% about image name
config.imgfold = fullfile(pwd,'demodata');
config.imgstr = 'A';
config.imgnum = 0:499;

%% set global var
global nodes mask filter_s rot_axis theta_delta
board = imread('board.jpg');
nodes = load('nodes.mat');nodes = nodes.nodes;
mask = load('mask.mat');mask = mask.mask;
filter_s = load('filter_s.mat');filter_s = filter_s.filter_s;
rot_axis = load('rot_axis.mat');rot_axis = rot_axis.rot_axis;
k = rot_axis.k;
b = rot_axis.b;
bottom = rot_axis.bottom;
theta_delta = 2*pi/length(config.imgnum);

%% 
if config.mask_dynamc_adj == true
    mask.dy_adj = true;
end

if config.write_into_txt == true
    filename = datestr(datetime);
    filename(filename==':') = '_';
    fp=fopen([filename,'.txt'],'a');
end

if config.realtime_disp == true
    figure
end

if config.read_img_prev == true
    % estimate memory for feasibility
    mem_img = numel(board)/size(board,3);
    meminfo = memory;
    mem_avi = meminfo.MemAvailableAllArrays;
    mem_need = length(config.imgnum) * mem_img;
    disp(['Need ',num2str(mem_need/2^30),...
        'GBytes to read all images in']);
    if mem_need + 2^30 > mem_avi
        error([...
            'Avilable memory is not enough to read all images in.',...
            ' Please change the option ''config.readallimgprev'' ',...
            'to ''false''']);
    end
    % init global big images mat
    global_img_mat('init',size(board,1),size(board,2),...
        length(config.imgnum),'uint8');
    % write into global big images mat
    for img_i = 1:length(config.imgnum)
        img_name = fullfile(config.imgfold,[config.imgstr,...
            num2str(config.imgnum(img_i)),'.jpg']);
        global_img_mat('set',img_i,img_name);
    end
end

%% main loop
disp('start processing')
if config.timing == true
    time_start = clock;
end

for img_i = 1:length(config.imgnum)
    % count
    if mod(img_i,50) == 0
        disp(num2str(img_i))
    end
    
    if config.read_img_prev == true
        img_name = img_i;
    else
        img_name = fullfile(config.imgfold,[config.imgstr,...
            num2str(config.imgnum(img_i)),'.jpg']);
    end
    
    % barycenter method to get centerline in images
    [m,n] = centerline_barycenter(img_name);
    % based on idx in images, calculate idx in chechkerboard
    [x,y] = real_position([m,n]);

    intersec_x = (k.*y+x-k*b)/(k^2+1);
    intersec_y = k.*intersec_x+b;
    distance = ((x-intersec_x).^2+(y-intersec_y).^2).^(1/2);
    
    coord_x = distance.*cos(theta_delta*img_i);
    coord_y = distance.*sin(theta_delta*img_i);
    coord_z = ((bottom(1)-intersec_x).^2+(bottom(2)-intersec_y).^2).^(1/2);
    
    
    if config.realtime_disp == true
        scatter3(coord_x,coord_y,coord_z,1,[1,0,0],'filled')
        set(gca,'DataAspectRatio',[20,20,20]),axis off,hold on
        pause(0.01)
    end

    if config.write_into_txt == true
        file_temp = [coord_x,coord_y,coord_z];
        for file_i = 1:size(file_temp,1)
            fprintf(fp,'%f %f %f\r\n',file_temp(file_i,:));
        end
    end

end

if config.timing == true
    total_time = etime(clock,time_start);
    aver_time = total_time/length(config.imgnum);
    disp(['total time: ',num2str(total_time),'s'])
    disp(['average time: ',num2str(aver_time),'s'])
end

if config.write_into_txt == true
    fclose(fp);
end

end