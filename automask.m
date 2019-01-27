function [bw] = automask(points_graph_c,board_ori)
board = rgb2gray(board_ori);
% figure,imshow(board)
[b_m,b_n] = size(board);

node_size = size(points_graph_c,2);
points_graph_c_ = points_graph_c(points_graph_c~=0);
points_graph_c_ = reshape(points_graph_c_,[],node_size);
[n_m,n_n] = size(points_graph_c_);
node_ld = points_graph_c_(n_m,1);
node_rd = points_graph_c_(n_m,n_n);
node_lu = points_graph_c_(1,1);
node_ru = points_graph_c_(1,n_n);

node_mat = complex2imgmat([node_ld;node_rd;node_ru;node_lu]);
y = node_mat(:,1)';
x = node_mat(:,2)';
bw = poly2mask(x,y,b_m,b_n);
% figure,imshow(bw)

figure,imshow(board,[])
hold on
for i = [node_ld,node_rd,node_lu,node_ru]
    rectangle('Position',[real(i)-26,-imag(i)-26,51,51],...
        'LineWidth', 3,'edgecolor','g','Curvature',[1 1])
end
x = [x,x(1)];
y = [y,y(1)];
line(x,y,'Color','red','LineWidth',2);
title('region of automask')
end
