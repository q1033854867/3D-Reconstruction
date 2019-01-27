function [nodes_neib] = find_node_neib(nodes_graph_c_r,nodes_graph_c_nz)
nodes_num = length(nodes_graph_c_r);
% node_neib = [upper lower left right upperleft loweright null null]
nodes_neib = zeros(nodes_num,8);
for node_i = 1:nodes_num
    cur_node = nodes_graph_c_r(node_i);
    [row,col] = find(nodes_graph_c_nz==cur_node);
    % upper
    try
        nodes_neib(node_i,1) = nodes_graph_c_nz(row-1,col);
    catch
        nodes_neib(node_i,1) = cur_node + 100i;
    end
    % lower
    try
        nodes_neib(node_i,2) = nodes_graph_c_nz(row+1,col);
    catch
        nodes_neib(node_i,2) = cur_node - 100i;
    end
    % left
    try
        nodes_neib(node_i,3) = nodes_graph_c_nz(row,col-1);
    catch
        nodes_neib(node_i,3) = cur_node - 100;
    end
    % right
    try
        nodes_neib(node_i,4) = nodes_graph_c_nz(row,col+1);
    catch
        nodes_neib(node_i,4) = cur_node + 100;
    end
    
    % upperleft
    try
        nodes_neib(node_i,5) = nodes_graph_c_nz(row-1,col-1);
    catch
        % node_neib(i,5) = cur_node + 100i;
    end
    % loweright
    try
        nodes_neib(node_i,6) = nodes_graph_c_nz(row+1,col+1);
    catch
        % node_neib(i,6) = cur_node - 100i;
    end
    
    % % left
    % try
    %     node_neib(i,7) = nodes_graph_c_nz(row,col-1);
    % catch
    %     node_neib(i,7) = cur_node - 100;
    % end
    % % right
    % try
    %     node_neib(i,8) = nodes_graph_c_nz(row,col+1);
    % catch
    %     node_neib(i,8) = cur_node + 100;
    % end
end
end

function [extend_neib] = node_graph_extend(nodes_graph_c,node,orient)
if nargin < 3
    % extend_neib = [upper lower left right]
    extend_neib = zeros(length(node),4);
    for node_i = 1:length(node)
        cur_node = node(node_i);
        [row,col] = find(nodes_graph_c==cur_node);
        % upper
        try
            extend_neib(node_i,1) = nodes_graph_c(row-1,col);
        catch
            extend_neib(node_i,1) = cur_node + 100i;
        end
        % lower
        try
            extend_neib(node_i,2) = nodes_graph_c(row+1,col);
        catch
            extend_neib(node_i,2) = cur_node - 100i;
        end
        % left
        try
            extend_neib(node_i,3) = nodes_graph_c(row,col-1);
        catch
            extend_neib(node_i,3) = cur_node - 100;
        end
        % right
        try
            extend_neib(node_i,4) = nodes_graph_c(row,col+1);
        catch
            extend_neib(node_i,4) = cur_node + 100;
        end
    end
else
    orient_num = length(orient);
    node_num = length(node);
    extend_neib = zeros(node_num,orient_num);
    for node_i = 1:node_num
        cur_node = node(node_i);
        [row,col] = find(nodes_graph_c==cur_node);
        for orient_i = 1:orient_num
            switch char(orient(orient_i))
                case 'left'
                    extend_neib(node_i,orient_i) = ...
                        nodes_graph_c(row,col-1);
                case 'upper'
                    extend_neib(node_i,orient_i) = ...
                        nodes_graph_c(row-1,col);
                case 'upperleft'
                    extend_neib(node_i,orient_i) = ...
                        nodes_graph_c(row-1,col-1);
                case 'loweright'
                    extend_neib(node_i,orient_i) = ...
                        nodes_graph_c(row+1,col+1);
            end
        end
    end
end
end