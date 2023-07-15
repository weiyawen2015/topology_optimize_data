clear all
close all
clc

% Img = imread('aaa.bmp');
% Imag = rgb2gray(Img);% 图像读取与灰度化
%  doubleGray = im2double(Imag);% 将图像数据uint8类型转换为double类型并进行归一化
%  A = doubleGray;
load S
x = 1:1:size(S,2);
for i = 1:size(S,2)
    CC = S{i};
    A  = dct2(CC);
    
   y(i,1) = A(1,1);
%      y(i,1) = max(max(A));
      yy(i,1) = min(min(A));
end

plot(x,y,'bo')
hold on
plot(x,yy,'r*')
hold off



% 
% 
% 
% 
% 
%  percent = 0.2;
%  
%  [nely nelx] = size(A);
%  
%  
% B = dct2(A);
%         
%  Bnew = B(1:floor(percent*nelx),1:floor(percent*nely));     
% huan= idct2(Bnew,[nely nelx] );
% 
%  
% imshow(log(abs(B)),[]),
% colormap(jet(64)),
% colorbar
%  title('频谱')
% 
%   figure(2)
% imshow(huan,'InitialMagnification','fit')
%   title('还原')
