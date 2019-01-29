function [m,n] = rotation_axis(img_addr,mask,varargin)
if nargin == 2
    % extract centerline method
    [m,n] = centerline_barycenter(img_addr,mask);
else
    % iteration method
    error('iteration method are under development')
end

end


