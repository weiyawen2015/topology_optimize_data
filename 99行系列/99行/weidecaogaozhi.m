clear all 
clc
close all
% global StrObj % Shared between myfun and myHessianFcn
% [StrObj] = StrSys();

x =-0.5:0.01:1.5;


alfa=3;

fxpar=[1e-9 1 2 3 ];

xpar=[0.5 1.5 2.5 ];


% flag==1
  tmpfx=(x-xpar(1));
 fx2=(1+tanh(alfa*tmpfx))/2 ;

   figure(2)
 plot(x,fx2,'b*')   
 
 
 
 
%  
%  
%    fx3= 0.5*(1+tanh(alfa*( x-0.9)));
%    
%     fx4 =  0.5*alfa*power(1-tanh(alfa*( x-0.9)),2);
%    
%  figure(2)
%  plot(x,fx2,'b*')   
%  
%         
%  figure(3)
%  plot(x,fx3,'b*')  
%  
%  figure(4)
%  plot(x,fx4,'b*')  
%  
% %% flag==2
% tmpfx=(x-xpar(1));
%  fx3= 1./(1+exp(-alfa*tmpfx));
%  
%   fx4 = (1./(1+exp(-alfa*(x-0.5))));
%  figure(3)
%  plot(x,fx4,'b*')     
%         
%         
% 
% pp=3;
% alfa=16;
% 
% fxpar=[1e-9 1 2 3 ];
% 
% xpar=[0.5 1.5 2.5 ];
% 
% for idx=1:length(xpar)
% 
%         filfun= (1+tanh(alfa*(x-xpar(idx))))/2 ;
% end