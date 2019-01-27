function [points_graph_c] = node_graph_baseright(points_raw,node_size,board)
% node number
num = size(points_raw);
points_raw = fliplr(points_raw);
%convert to complex
points_raw = points_raw(:,1)-1i.*points_raw(:,2);

img_rightbottom = size(board,2)-1i.*size(board,1);
distance = abs(points_raw-img_rightbottom);
[~,idx] = min(distance);
cal = points_raw(idx);

points_graph_c = zeros(node_size,node_size);
points_graph_c(node_size,node_size) = cal;

% figure,imshow(board)
% hold on

% ×ó ×óÉÏ ÉÏ
base_m = [0,-1,-1];
base_n = [-1,-1,0];

points_raw_ = points_raw;

for row = node_size:-1:1
    for col = node_size:-1:1
        cur_p = points_graph_c(row,col);
        
%         rectangle('Position',[real(cur_p)-16,-imag(cur_p)-16,31,31],...
%             'edgecolor','g','Curvature',[1 1])
        
        points_raw_(points_raw_==cur_p) = [];
        
        delta = points_raw_ - cur_p;
        delta_ang = angle(delta);
        delta_ang(delta_ang<0) = 2*pi + delta_ang(delta_ang<0);
        
        % angle constraint
        % subject to left nodes
        ang_range = pi.*[7/6,3/4];
        ang_constraint_idx = ...
            delta_ang<ang_range(1) & delta_ang>ang_range(2);
        % num
        % disp(sum(ang_constraint_idx))
        if sum(ang_constraint_idx) == 0
            continue
        end
        index_constraint = points_raw_.*ang_constraint_idx;
        index_constraint(index_constraint==0) = [];
        delta_constraint = index_constraint - cur_p;
        % distance ascent
        [~,d_ix] = sort(abs(delta_constraint));
        points_graph_c(row+base_m(1),col+base_n(1)) = ...
            index_constraint(d_ix(1));
        
%         pause(0.1)
    end
    
    % row_start
    cur_p = points_graph_c(row,node_size);
    delta = points_raw_ - cur_p;
    delta_ang = angle(delta);
    delta_ang(delta_ang<0) = 2*pi + delta_ang(delta_ang<0);
    % constraint to up nodes
    ang_range = pi.*[3/4,1/4];
    ang_constraint_idx = ...
        delta_ang<ang_range(1) & delta_ang>ang_range(2);
    % num
    % disp(sum(ang_constraint_idx))
    if sum(ang_constraint_idx) == 0
        break
    end
    index_constraint = points_raw_.*ang_constraint_idx;
    index_constraint(index_constraint==0) = [];
    delta_constraint = index_constraint - cur_p;
    % distance ascent
    [~,d_ix] = sort(abs(delta_constraint));
    points_graph_c(row+base_m(3),17+base_n(3)) = ...
        index_constraint(d_ix(1));
end

end

