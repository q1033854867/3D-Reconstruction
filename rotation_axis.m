function [m,n] = rotation_axis(img_addr,mask)
% extract centerline method
[m,n] = centerline_barycenter(img_addr,mask);
% iteration method
end

