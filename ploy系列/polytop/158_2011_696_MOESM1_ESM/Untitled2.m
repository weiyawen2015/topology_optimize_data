close all
clear all
clc
%%
step1 = 2;
load resultmesh
result = resultmesh;
for i = 1:step1:size(result,2)
  img = result{i};
 img = imresize(img,1);
%-�ڴ�-------------------------------
[I,map]=rgb2ind(img,256); 
%--------------------------------   
    if i==1
        imwrite(I,map,'Gif_polymesh.gif','gif', 'Loopcount',inf,'DelayTime',0.000001);%��һ�α��봴����
    elseif i==size(result,2)
        imwrite(I,map,'Gif_polymesh.gif','gif','WriteMode','append','DelayTime',0.000001);
    else
        imwrite(I,map,'Gif_polymesh.gif','gif','WriteMode','append','DelayTime',0.000001);
    end
    i
end