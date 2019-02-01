% clear all, close all
load('img_crop.mat')
[col_max_val,col_max_idx] = max(img_crop,[],1);
figure,imshow(img_crop,[])
hold on
for j = 1:size(img_crop,2)
    scatter(j,col_max_idx(j),1,'filled');
end