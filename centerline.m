function [m,n] = centerline(add)
% ����
% ����mask��ĻҶ�ͼ
% addΪ��ַ
% ���
% mΪ�������ĵ�������
% nΪ�������ĵ�������

h1 = fspecial('average',[5,15]);
%�˲���size
h2_height = 15;
h2_width = 31;
h2 = -1/(h2_height*h2_width).*ones(h2_height,h2_width);
h2((h2_height+1)/2,(h2_width+1)/2) = 1;

img = imread(add);
img = im2double(img);

img = imfilter(img,h1,'conv','symmetric','same');
img = imfilter(img,h2,'conv','symmetric','same');

img(find(img < 0.1*max(img(:)))) = 0;

[max_value index] = max(img');
img = 0.*img;

for j = 1:1200
    if max_value(j) > 0
        img(j,index(j)) = 1;
    end
end
[m,n] = find(img == 1);

