clear all
close all
clc

% Img = imread('aaa.bmp');
% Imag = rgb2gray(Img);% ͼ���ȡ��ҶȻ�
%  doubleGray = im2double(Imag);% ��ͼ������uint8����ת��Ϊdouble���Ͳ����й�һ��
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
%  title('Ƶ��')
% 
%   figure(2)
% imshow(huan,'InitialMagnification','fit')
%   title('��ԭ')
