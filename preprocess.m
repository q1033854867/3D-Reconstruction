function preprocess(varargin)
% usage:

if nargin == 0
    board_path = fullfile(pwd,'board.jpg');
    rot_axis_path = fullfile(pwd,'rot_axis.jpg');
    
    nodes = get_nodes(board_path);
    save('nodes.mat','nodes')
    
    mask = get_mask(nodes,board_path,rot_axis_path);
    save('mask.mat','mask')
    
    rot_axis = get_rot_axis(nodes,mask,rot_axis_path);
    save('rot_axis.mat','rot_axis')
    
    filter_s = get_filter_s();
    save('filter_s.mat','filter_s')
    
else

    if sum(strcmp(varargin,'board_path')) == 1
        board_path_idx = find(strcmp(varargin,'board_path'));
        board_path = varargin{board_path_idx+1};
    else
        board_path = fullfile(pwd,'board.jpg');
    end

    if sum(strcmp(varargin,'rot_axis_path')) == 1
        rot_axis_path_idx = find(strcmp(varargin,'rot_axis_path'));
        rot_axis_path = varargin{rot_axis_path_idx+1};
    else
        rot_axis_path = fullfile(pwd,'rot_axis.jpg');
    end

    if sum(strcmp(varargin,'filter_s')) == 1
        filter_s_idx = find(strcmp(varargin,'filter_s'));
        filter_s_var1 = varargin{filter_s_idx+1};
        filter_s_var2 = varargin{filter_s_idx+2};
        filter_s = get_filter_s(filter_s_var1,filter_s_var2);
        save('filter_s.mat','filter_s')
    end

end
end