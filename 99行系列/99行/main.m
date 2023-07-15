clear all
close all
clc
nelx=60;   % x轴方向的单元数
nely=20;   % y轴方向单元数
volfrac=0.5;  %体积比
% penal=3;      %材料插值的惩罚因子
rmin=2.4;       %敏度过滤的半径
alfa=2.5;
xfar = 1;




xshi =0:0.01:1;
% flag==1
  tmpfx=(xshi-xfar);
 fx=(1+tanh(alfa*tmpfx));
 fx2 = xshi.^3;
  figure(1)
  plot(xshi,fx,'b*')
  hold on
   plot(xshi,fx2,'ro') 
  hold off






% INITIALIZE
x(1:nely,1:nelx) = volfrac; % x是设计变量 赋初始值
loop = 0; 
change = 1.;
% START ITERATION
 while change > 0.01  

  loop = loop + 1;
  xold = x; % 存储前一次变量
  
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,alfa,xfar);      %进行有限元分析      
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;  %单元刚度矩阵
  
  c = 0.; %用来存放目标函数的变量.这里目标函数是刚度最大，也就是柔
                                     %度最小  
   % 单元挨个扫描                                  
  for ely = 1:nely  
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);% 8x1的矩阵单元的位移向量
%       c = c + x(ely,elx)^penal*Ue'*KE*Ue;       %计算目标函数的值（即柔度是多少）
% 
        c = c + (1+tanh(alfa*( x(ely,elx)-xfar)))*Ue'*KE*Ue;       %计算目标函数的值（即柔度是多少
        dc(ely,elx) = -alfa*(1-tanh(alfa*( x(ely,elx)-xfar))^2)*Ue'*KE*Ue; % 灵敏度分析的结果 为设计变量更新做准备
   
        
%         c = c + (1/(1+exp(-alfa*(x(ely,elx)-0.5))))*Ue'*KE*Ue;       %计算目标函数的值（即柔度是多少
%         dc(ely,elx) = -alfa*exp(-alfa*(x(ely,elx)-0.5))/(1+exp(-alfa*(x(ely,elx)-0.5)))^2*Ue'*KE*Ue; % 灵敏度分析的结果 为设计变量更新做准备     
    
    end
  end
  
  
% FILTERING OF SENSITIVITIES
   [dc]   = check(nelx,nely,rmin,x,dc);          %灵敏度过滤，为了边界光顺一点
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dc);       %优化准则法更新设计变量 这里是xnew
  save x x
% PRINT RESULTS
  change = max(max(abs(x-xold))); %新的设计变量与旧的设计变量的差的绝对值得最大值
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])          %屏幕上显示迭代信息
% PLOT DENSITIES  
  figure(2)
  colormap(gray); % x是矩阵 不是索引图像 加"-"是为了黑白颠倒一下
  imagesc(-x);  %colormap(gray);imagesc(-x);配合使用 详细解释看https://blog.csdn.net/zy122121cs/article/details/49761307
  axis equal; 
  axis tight; % 显示紧凑一些
  axis off;
  pause(1e-7); 
%    title([' 第loop次循环' num2str(loop)]) 
  %优化结果的图形显示（个人认为这种图形显示方法很不好，太简单了。比较方便的图形显示应该是：
% 每一次迭代同时显示优化结果、目标函数曲线，然后自动保存每一次的结果）

end
