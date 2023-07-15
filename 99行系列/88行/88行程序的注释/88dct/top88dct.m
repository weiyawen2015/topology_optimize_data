%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%


clear 
clc
nelx=120;   % x轴方向的单元数
nely=40;   % y轴方向单元数

vvnelx=12;
uunely=5;
volfrac=0.5;  %体积比
penal=3;      %材料插值的惩罚因子
alfa = 20;
%  function top88(nelx,nely,volfrac,penal,rmin,ft)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);%单元刚度阵
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);%(nely+1)*(nelx+1)个节点
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);%取并
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);%取减
%% 灵敏度准备
  alfa_u = sqrt(2/nely)*ones(uunely,1);
   alfa_u(1,1) = sqrt(1/nely);
    alfa_v = sqrt(2/nelx)*ones(vvnelx,1);
   alfa_v(1,1) = sqrt(1/nelx);
   
   for u = 1:uunely
      for v = 1:vvnelx
           dxT_x = zeros(nely,nelx);
          for i = 1:nely
               for j = 1:nelx
                   dxT_x(i,j) = alfa_u(u,1)*alfa_v(v,1)*cos(0.5*pi*(2*(i-1)+1)*(u-1)/nely)*cos(0.5*pi*(2*(j-1)+1)*(v-1)/nelx);
               end
          end
          S{u,v} =  dxT_x;
      end
   end

%% INITIALIZE ITERATION iteration循环 迭代
 x_ori = repmat(volfrac,nely,nelx);%给个初始值  随便给
[xxinit]=dct2(x_ori);
x = xxinit(1:uunely,1:vvnelx);%repmat(1,uunely,vvnelx)
xTilde = idct2(reshape(x,uunely,vvnelx),nely,nelx);
xPhys = (1+tanh(alfa*(xTilde-0.5)))/2;
loop = 0;
change = 1;

%% START ITERATION
% while change > 0.01
while loop < 10
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);%解有限元方程
  %% objective function and sensitivity analysis目标函数 和灵敏度分析
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)); % 目标函数
  %---------------------------------------------------------
  %灵敏度分析
   dc_rou = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
   drou_xT = 0.5*alfa*(1-tanh(alfa*(xTilde-0.5)).^2);
   dv = ones(uunely,vvnelx);%V对设计变量的灵敏度一直是1
   dc = zeros(uunely,vvnelx);
    for u = 1:uunely
       for v = 1:vvnelx
            dc(u,v) = sum(sum(dc_rou.*drou_xT.*S{u,v}));
%             dv(u,v) = sum(sum(dv_rou.*drou_xT.*S{u,v}));
       end
    end
  %---------------------------------------------------------
  %% OC准则 OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES criterion标准variables变量 
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
     x = xnew;
     xTilde=   idct2(reshape(x,uunely,vvnelx),nely,nelx);
     xPhys = (1+tanh(alfa*(xTilde-0.5)))/2;
    % 二分法
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
