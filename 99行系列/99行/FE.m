function [U]=FE(nelx,nely,x,alfa,xfar)
[KE] = lk; %单元刚度矩阵
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
    %总体刚度矩阵的稀疏矩阵 注意节点数和单元数的关系
F = sparse(2*(nely+1)*(nelx+1),1); %载荷向量 有这么多节点(nely+1)*(nelx+1)
U = zeros(2*(nely+1)*(nelx+1),1); %位移向量 有这么多节点(nely+1)*(nelx+1)
% 总刚组装
for elx = 1:nelx   % 第elx行的单元
  for ely = 1:nely  %第ely列的单元
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];%这里的Y轴是反向的，但是不影响最后的结果，详情请见二楼TYNGOD这位高手的解释，感谢TYNGOD。
%     K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;  %将单元刚度矩阵组装成总的刚度矩阵
       K(edof,edof) = K(edof,edof) +  (1+tanh(alfa*( x(ely,elx)-xfar)))*KE;  %将单元刚度矩阵组装成总的刚度矩阵
%             K(edof,edof) = K(edof,edof) + (1/(1+exp(-alfa*(x(ely,elx)-0.5))))*KE;  %将单元刚度矩阵组装成总的刚度矩阵
  
  
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F(2,1) = -1; % 其他元素都是0

fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]); %union取并运算固定结点
alldofs     = [1:2*(nely+1)*(nelx+1)];                       %所有结点 
freedofs    = setdiff(alldofs,fixeddofs);  % setdiff返回A中有，B中没有的值自由节点
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:); %求解   
U(fixeddofs,:)= 0;












