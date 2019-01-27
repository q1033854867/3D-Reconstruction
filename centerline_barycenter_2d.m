function [m,n] = centerline_barycenter_2d(img_addr,mask,h1_size,h2_size)
% 输入
% addr为原图像地址
% mask_v有效点模板
% 输出
% m为光条中心点行坐标
% n为光条中心点列坐标
% t1 = clock;
if nargin < 3
    h1_size = 9;
    h2_size = 17;
end

img = rgb2gray(imread(img_addr));
img = im2double(img);
img = img.*mask;
[img_m,img_n] = size(img);

col = find(sum(mask==1,1));
row = find(sum(mask==1,2));
c_min = min(col);
c_max = max(col);
r_min = min(row);
r_max = max(row);
img_crop = img(r_min:r_max,c_min:c_max);

h1 = fspecial('average',[h1_size,h1_size]);
%滤波器size
h2_height = h2_size;
h2_width = h2_size;
h2 = -1/(h2_height*h2_width).*ones(h2_height,h2_width);
h2((h2_height+1)/2,(h2_width+1)/2) = 1;

img_crop = imfilter(img_crop,h1,'conv','symmetric','same');
img_crop = imfilter(img_crop,h2,'conv','symmetric','same');

% t2 = clock;
% etime(t2,t1)
% t1 = clock;

img = img*0;
img(r_min:r_max,c_min:c_max) = img_crop;

img(img < 0.1*max(img(:))) = -1;

m = [];
n = [];

[max_value index] = max(img');

for j = r_min:r_max
    if max_value(j) > 0
        dev = 1;
        sum_v = max_value(j);
        w_sum = sum_v*index(j);
        th = 0.1*max_value(j);
        while img(j,index(j) + dev) > th
            sum_v = sum_v + img(j,index(j) + dev);
            w_sum = w_sum + (index(j) + dev)*img(j,index(j) + dev);
            dev = dev + 1;
        end
        dev = 1;
        while img(index(j) - dev) > th
            sum_v = sum_v + img(index(j) - dev);
            w_sum = w_sum + (index(j) - dev)*img(j,index(j) - dev);
            dev = dev + 1;
        end
        w_cen = w_sum/sum_v;
        m = [m;j];
        n = [n;w_cen];
    end
end

% m = [];
% n = [];

img = imrotate(img,-90);
[max_value index] = max(img');
for j = c_min:c_max
    if max_value(j) > 0
        dev = 1;
        sum_v = max_value(j);
        w_sum = sum_v*index(j);
        th = 0.1*max_value(j);
        while img(j,index(j) + dev) > th
            sum_v = sum_v + img(j,index(j) + dev);
            w_sum = w_sum + (index(j) + dev)*img(j,index(j) + dev);
            dev = dev + 1;
        end
        dev = 1;
        while img(index(j) - dev) > th
            sum_v = sum_v + img(index(j) - dev);
            w_sum = w_sum + (index(j) - dev)*img(j,index(j) - dev);
            dev = dev + 1;
        end
        w_cen = w_sum/sum_v;
%         m = [m;j];
%         n = [n;w_cen];
        m = [m;img_m-w_cen+1];
        n = [n;j];
    end
end
% t2 = clock;
% etime(t2,t1)

end
