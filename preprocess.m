clear all,close all
board = imread(fullfile(pwd,'board.jpg'));
board_ori = board;
board = im2double(rgb2gray(board));

board = board.^2;
% imhist(board)
% pause

board = double(im2bw(board,0.4));
% figure,imshow(board)
% pause

% corr base
base_size = 60;
base = zeros(base_size)-1;
base(1:round(base_size/2),1:round(base_size/2)) = 1;
base(round(base_size/2)+1:end,round(base_size/2)+1:end) = 1;

nodes = imfilter(board,base,'conv','symmetric','same');

th = 0.6;
cor_min = min(nodes(:));
cor_max = max(nodes(:));
nodes(nodes<th*cor_max & nodes>th*cor_min) = 0;
nodes(nodes~=0) = 1;
% corr result
% figure,imshow(points,[])
% pause

single_points = nodes*0;
L = bwlabel(nodes);
count = tabulate(L(:));
count = count(2:end,:);

nodes_raw = zeros(size(count));
nodes_raw = nodes_raw(:,1:2);
for i = 1:size(count,1)
    if count(i,2) <= 10
        nodes(find(L==i)) = 0;
    else
        [m,n] = find(L==i);
        x = int32((max(n)+min(n))/2);
        y = int32((max(m)+min(m))/2);
        nodes_raw(i,:) = [y x];
%         single_points(y,x) = 1;
    end
end
nodes_raw(all(nodes_raw==0,2),:) = [];

figure,imshow(board_ori)
for i = 1:size(nodes_raw,1)
    rectangle('Position',[nodes_raw(i,2)-16,nodes_raw(i,1)-16,31,31],...
        'edgecolor','g','Curvature',[1 1])
end
title('checkerboard nodes')
pause(0.1)

% figure,imshow(single_points)
% hold on
% scatter(index(:,2),index(:,1),5,[1,1,1],'s','filled')
% for i = 1:size(index,1)
%     rectangle('Position',[index(i,2)-16,index(i,1)-16,31,31],...
%         'edgecolor','g','Curvature',[1 1])
% end
% pause(0.1)

if ~exist(fullfile(pwd,'nodes_raw.mat'))
    save('nodes_raw.mat','nodes_raw')
end

node_size = 17;
nodes_graph_c = node_graph_baseright(nodes_raw,node_size,board_ori);
if ~exist(fullfile(pwd,'nodes_graph_c.mat'))
    save('nodes_graph_c.mat','nodes_graph_c')
end

nodes_graph_c_r = nodes_graph_c(nodes_graph_c~=0);
nodes_graph_c_nz = reshape(nodes_graph_c_r,[],node_size);
nodes_neib = find_node_neib(nodes_graph_c_r,nodes_graph_c_nz);
if ~exist(fullfile(pwd,'nodes_graph_c_r.mat'))
    save('nodes_graph_c_r.mat','nodes_graph_c_r')
end
if ~exist(fullfile(pwd,'nodes_graph_c_nz.mat'))
    save('nodes_graph_c_nz.mat','nodes_graph_c_nz')
end
if ~exist(fullfile(pwd,'nodes_neib.mat'))
    save('nodes_neib.mat','nodes_neib')
end

mask = automask(nodes_graph_c,board_ori);
if ~exist(fullfile(pwd,'mask.mat'))
    save('mask.mat','mask')
end

mask_mini = automask_mini(nodes_graph_c,board_ori);
if ~exist(fullfile(pwd,'mask_mini.mat'))
    save('mask_mini.mat','mask_mini')
end

rot_axis_filename = fullfile(pwd,'rot_axis.jpg');
[m,n] = rotation_axis(rot_axis_filename,mask_mini);
rot_axis.img_idx = [m,n];
rot_axis_idx = real_position([m,n],nodes_graph_c_r,...
    nodes_graph_c_nz,nodes_neib);
rot_axis.real_idx = rot_axis;
p = polyfit(rot_axis_idx(:,1),rot_axis_idx(:,2),1);
rot_axis.k = p(1);
rot_axis.b = p(2);
% k = p(1);
% b = p(2);
rot_axis.bottom = rot_axis_idx(1,:);
mask_mini_rotaxis = automask_mini_rotaxis(nodes_graph_c,board_ori,[m,n]);
if ~exist(fullfile(pwd,'mask_mini_rotaxis.mat'))
    save('mask_mini_rotaxis.mat','mask_mini_rotaxis')
end
if ~exist(fullfile(pwd,'rot_axis.mat'))
    save('rot_axis.mat','rot_axis')
end