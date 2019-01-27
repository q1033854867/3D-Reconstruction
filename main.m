clear all,close all

config.writeintotxt = false;
config.realtimedisp = true;
config.readallimgprev = true
config.timing = true;
config.imgfold = fullfile(pwd,'demodata');
config.imgstr = 'A';
config.imgnum = 0:499;

load('rot_axis.mat')
load('mask_mini_rotaxis.mat')
load('nodes_graph_c_r.mat')
load('nodes_graph_c_nz.mat')
load('nodes_neib.mat')
board = imread('board.jpg');

k = rot_axis.k;
b = rot_axis.b;
bottom = rot_axis.bottom;

theta_delta = 2*pi/length(config.imgnum);

if config.writeintotxt == true
    filename = datestr(datetime);
    filename(filename==':') = '_';
    fp=fopen([filename,'.txt'],'a');
end
if config.realtimedisp == true
    figure
end
if config.timing == true
    time_start = clock;
end

ttt = zeros(1,3);
global tttt
tttt = zeros(1,4);

for img_i = config.imgnum
    
    if mod(img_i,20) == 0
        disp(img_i)
    end
    
    img_name = fullfile(config.imgfold,...
        [config.imgstr,num2str(img_i),'.jpg']);
    
    t1 = clock;
    [m,n] = centerline_barycenter(img_name,mask_mini_rotaxis,...
        [17,17],[31,31]);
    ttt(1) = ttt(1) + etime(clock,t1);
    
    t1 = clock;
    laser_idx = real_position([m,n],nodes_graph_c_r,...
        nodes_graph_c_nz,nodes_neib);
    ttt(2) = ttt(2) + etime(clock,t1);
    
    
    t1 = clock;
    x = laser_idx(:,1);
    y = laser_idx(:,2);
    
    % 此处的k与b为旋转中轴的
    intersec_x = (k.*y+x-k*b)/(k^2+1);
    intersec_y = k.*intersec_x+b;
    distance = ((x-intersec_x).^2+(y-intersec_y).^2).^(1/2);
    
    coord_x = distance.*cos(theta_delta*img_i);
    coord_y = distance.*sin(theta_delta*img_i);
    coord_z = ((bottom(1)-intersec_x).^2+(bottom(2)-intersec_y).^2).^(1/2);
    
    ttt(3) = ttt(3) + etime(clock,t1);
    
    if config.realtimedisp == true
        scatter3(coord_x,coord_y,coord_z,1,[1,0,0],'filled')
        set(gca,'DataAspectRatio',[20,20,20])
        axis off
        hold on
        pause(0.01)
    end
    if config.writeintotxt == true
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
if config.writeintotxt == true
    fclose(fp);
end