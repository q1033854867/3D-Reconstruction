function [m,n] = centerline_barycenter(img_addr,mask,h1_size,h2_size)
% 输入
% addr为原图像地址
% mask_v有效点模板
% 输出
% m为光条中心点行坐标
% n为光条中心点列坐标
% t1 = clock;
if nargin < 3
    h1_size = [11,11];
    h2_size = [21,21];
end

% global tttt
% t1 = clock;

img = imread(img_addr);
img = rgb2gray(img);
img = im2double(img);
img = img.*mask;

% tttt(1) = tttt(1) + etime(clock,t1);
% t1 = clock;

col = find(sum(mask==1,1));
row = find(sum(mask==1,2));
c_min = min(col);
c_max = max(col);
r_min = min(row);
r_max = max(row);
img_crop = img(r_min:r_max,c_min:c_max);

% tttt(2) = tttt(2) + etime(clock,t1);
% t1 = clock;

h1 = fspecial('average',h1_size);
h2 = -1/(h2_size(1)*h2_size(2)).*ones(h2_size(1),h2_size(2));
h2((h2_size(1)+1)/2,(h2_size(2)+1)/2) = 1;
h = conv2(h1,h2);
% img_crop = imfilter(img_crop,h,'conv','symmetric','same');
img_crop = conv2(img_crop,h,'same');
% img_crop = imfilter(img_crop,h1,'conv','symmetric','same');
% img_crop = imfilter(img_crop,h2,'conv','symmetric','same');

% tttt(3) = tttt(3) + etime(clock,t1);
% t1 = clock;

% 移到前面
img_crop(img_crop < 0.1*max(img_crop(:))) = -1;
% search_points(img_crop);
[~,m,n] = search_points(img_crop);
m = m + r_min;
n = n + c_min;

% tttt(4) = tttt(4) + etime(clock,t1);


end
