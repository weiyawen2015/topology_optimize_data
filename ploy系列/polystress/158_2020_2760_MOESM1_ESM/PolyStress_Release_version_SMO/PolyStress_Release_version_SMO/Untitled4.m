close all
clear all
clc
%%
step1 = 1;
load result
for i = 1:step1:size(result,2)
  img = result{i};
 img = imresize(img,0.5);
%-内存-------------------------------
[I,map]=rgb2ind(img,256); 

%--------------------------------   
    if i==1
        imwrite(I,map,'Gif_polystress.gif','gif', 'Loopcount',inf,'DelayTime',0.000001);%第一次必须创建！
    elseif i==size(result,2)
        imwrite(I,map,'Gif_polystress.gif','gif','WriteMode','append','DelayTime',0.000001);
    else
        imwrite(I,map,'Gif_polystress.gif','gif','WriteMode','append','DelayTime',0.000001);
    end
    i
end