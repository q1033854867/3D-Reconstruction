function [real_idx] = real_position(points,nodes_graph_c_r,...
    nodes_graph_c_nz,nodes_neib,board)
% points: laserline points
% nodes_graph_c_r: checkerboard complex nodes(reshape to 1d)
% nodes_graph_c_nz: checkerboard complex nodes(without zero points)
% board: checkerboard image

% plot scatter
if nargin == 5
    node = complex2imgmat(nodes_graph_c_r);
    figure,imshow(board)
    hold on
    scatter(node(:,2),node(:,1),'filled')
    scatter(points(:,2),points(:,1),'filled')
    pause(0.01)
end

nodes_num = size(nodes_graph_c_r,1);
% [n_m,n_n] = size(nodes_graph_c_nz);

points_c = imgmat2complex(points);
points_num = size(points,1);

delta = points_c*ones(1,nodes_num) - ones(points_num,1)*nodes_graph_c_r.';
distance = abs(delta);
[~,n_nodes_idx] = min(distance,[],2);
n_nodes_idx_uni = unique(n_nodes_idx);

% nearest node
n_nodes = nodes_graph_c_r(n_nodes_idx);
n_nodes_uni = nodes_graph_c_r(n_nodes_idx_uni);
% 8 neighbors of nearest node
n_nodes_neib8 = nodes_neib(n_nodes_idx,:);
n_nodes_neib8_uni = nodes_neib(n_nodes_idx_uni,:);
n_nodes_neib4_uni = n_nodes_neib8_uni(:,1:4);

n_nodes_neib4_uni_angl = angle(n_nodes_neib4_uni - n_nodes_uni);
n_nodes_neib4_uni_angl_rot = ...
    rot_angle(n_nodes_neib4_uni_angl,n_nodes_neib4_uni_angl(:,4));

n_nodes_bw = (n_nodes*ones(1,length(n_nodes_uni)) == n_nodes_uni.');

n_nodes_neib4_angl = n_nodes_bw * n_nodes_neib4_uni_angl;
n_nodes_neib4_angl_rot = n_nodes_bw * n_nodes_neib4_uni_angl_rot;

n_nodes_angl = angle(points_c - n_nodes);
n_nodes_angl_rot = rot_angle(n_nodes_angl,n_nodes_neib4_angl(:,4));


n_nodes_orient = false(points_num,4);
% 最近节点在左下方
n_nodes_orient(:,1) = n_nodes_angl_rot<n_nodes_neib4_angl_rot(:,1);
% 最近节点在右下方
n_nodes_orient(:,2) = n_nodes_angl_rot<n_nodes_neib4_angl_rot(:,3)&...
    n_nodes_angl_rot>=n_nodes_neib4_angl_rot(:,1);
% 最近节点在右上方
n_nodes_orient(:,3) = n_nodes_angl_rot<n_nodes_neib4_angl_rot(:,2)&...
    n_nodes_angl_rot>=n_nodes_neib4_angl_rot(:,3);
% 最近节点在左上方
n_nodes_orient(:,4) = 1 - sum(n_nodes_orient(:,1:3),2);


% loweright node in neighbor rectangle
p_rec_neib = zeros(points_num,1);
% 最近节点在右下方 直接赋值
p_rec_neib(n_nodes_orient(:,2)) = n_nodes(n_nodes_orient(:,2));
% 最近节点在左下方 
p_rec_neib(n_nodes_orient(:,1)) = n_nodes_neib8(n_nodes_orient(:,1),4);
% 最近节点在右上方
p_rec_neib(n_nodes_orient(:,3)) = n_nodes_neib8(n_nodes_orient(:,3),2);
% 最近节点在左上方
p_rec_neib(n_nodes_orient(:,4)) = n_nodes_neib8(n_nodes_orient(:,4),6);

% 对于右下角点 顺序为自己 上 左上 左
p_rec_neib_uni = unique(p_rec_neib(:,1));

p_lr_idx = (1:nodes_num)*...
    (nodes_graph_c_r*ones(1,length(p_rec_neib_uni))==p_rec_neib_uni.');

% node_neib = [upper lower left right upperleft loweright null null]
orient_bw = logical([1,0,1,0,1,0,0,0]);
% 现在的顺序为 右下 右上 左下 左上
p_rec_neib_uni = [p_rec_neib_uni,nodes_neib(p_lr_idx,orient_bw)];
% change order
p_rec_neib_uni(:,[4,3]) = p_rec_neib_uni(:,[3,4]);

real_idx = zeros(points_num,2);
% loweright upperight upperleft lowerleft
A = [0,0;0,1;-1,1;-1,0];
for tf_i = 1:size(p_rec_neib_uni,1)
    B = complex2imgmat(p_rec_neib_uni(tf_i,:).');
    % TForm = cp2tform(B,A,'projective');
    tform = fitgeotrans(B,A,'projective');
    
    rl_bw = (p_rec_neib==p_rec_neib_uni(tf_i,1));
    
    cur_points = complex2imgmat(points_c(rl_bw));
    
    d_xy = transformPointsForward(tform,cur_points);
    % [d_x,d_y] = tformfwd(tform,cur_points);
    
    [node_y,node_x] = find(nodes_graph_c_nz==p_rec_neib_uni(tf_i,1));
    real_idx(rl_bw,:) = [node_x,node_y] + d_xy*diag([1,-1]);
    % real_idx(rl_bw,:) = [node_x,node_y] + [d_x,-d_y];
end

% projection transform result
if nargin == 5
    figure
    scatter(real_idx(:,1),real_idx(:,2),1,[1,0,0],'filled')
    axis on
    axis([0,17,0,17])
end
end

function [result] = rot_angle(input,base)
%
input(input<0) = 2*pi + input(input<0);
result = mod(input - base,2*pi);
end