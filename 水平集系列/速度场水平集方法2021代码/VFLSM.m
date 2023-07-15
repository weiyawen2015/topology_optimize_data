function VFLSM(nx,ny,vf,xnum) % A Matlab (2010a or later) code for VFLS
%% -------------------FE mesh & Level Set (LS) grid-----------------
nod_n = [ny+1,nx+1]; FEnod_glob = reshape(1:prod(nod_n),nod_n); % Node
FEnod_loc = reshape(FEnod_glob(1:ny,1:nx),1,nx*ny); % Local node index
FEnod_loc = repmat(FEnod_loc,4,1)+repmat([0;1;ny+2;ny+1],1,nx*ny);
LSgrid = nod_n+1; [LSy,LSx] = ndgrid(-0.5:ny+0.5,-0.5:nx+0.5); % LS grid
[fx,fy,fixed] = deal(nx/2+1,ny+1,[1,2,2*(ny+1)*nx+2]); % Boundary condition
f = sparse(2*(ny+1)*(fx-1)+2*fy,1,-1,2*(nx+1)*(ny+1),1); % Load
freedof = 1:2*prod(nod_n); freedof(fixed) = []; % Free node
FE_dof = 2*FEnod_loc(repmat(1:4,2,1),:)-repmat([1;0],4,nx*ny); % DOF index
[Ipp,Jpp] = deal(FE_dof(:,repmat(1:nx*ny,8,1)),FE_dof(repmat(1:8,8,1),:));
half = (Ipp(:)<=Jpp(:)); % Preparation for assembling stiffness matrix
%% --------------------------Gauss points---------------------------
gs = [(2*70^0.5)/63+5/9,5/9-(2*70^0.5)/63].^0.5; gs = [-gs,0,gs];
gw = 2./(1-gs.^2)./((315*gs.^4)/8-(105*gs.^2)/4+15/8).^2;
[gy,gx] = ndgrid(gs,gs); [gwy,gwx] = ndgrid(gw,gw);
gweight = gwx(:).*gwy(:); gnum = length(gweight); % Weighting & number
Dmat = (1/(1-0.3^2))*[1  0.3  0; 0.3  1  0; 0  0  (1-0.3)/2];
Bgs = @(i)0.5*[-1+gy(i),0 ,-gy(i)-1 ,0,1+gy(i),0 ,1-gy(i),0;
    0,-1+gx(i),0,1-gx(i),0,gx(i)+1,0,-gx(i)-1;
    -1+gx(i),-1+gy(i),1-gx(i),-gy(i)-1,gx(i)+1,1+gy(i),-gx(i)-1,1-gy(i)];
BDB = @(i)reshape(Bgs(i)'*Dmat*Bgs(i),64,1)*gweight(i)/4;
Gauss_k = cell2mat(arrayfun(@(i)BDB(i),1:gnum,'UniformOutput',false));
%% -------------------------Initialization--------------------------
[ynum,cr] = deal(max(floor(xnum/nx*ny),1),nx/xnum/4);
[cy,cx] = ndgrid((0:2:ynum*2)/(2*ynum/ny),(1:2:xnum*2+1)/(2*xnum/nx));
[cx,cy] = deal([cx(:);cx(:)-nx/xnum/2],[cy(:);cy(:)+ny/ynum/2]);
cNAN = max(abs([cx-fx,cy-fy])-1,[],2)<=cr;  [cx(cNAN),cy(cNAN)] = deal(inf);
[cx,xt] = ndgrid(cx,LSx);  [cy,yt] = ndgrid(cy,LSy);
Phi = sqrt((xt-cx).^2+(yt-cy).^2)-cr; % Initial LS function
Phi = reshape(min([Phi;[LSx(:),nx-LSx(:),LSy(:),ny-LSy(:)]'+0.5]),LSgrid);
Iter = 0;  obj_old = inf(4,1);  obj = 0;  con = 1;
[Vn,Vn_old1,Vn_old2,dobj,dcon] = deal(zeros(LSgrid)); % Des.var.& sensi.
[Vn_max,upp,low] = deal(0.5,0.5,-0.5); % Des. var. bounds & asymptote in MMA
cc = @(obj,con)obj/max((con+vf),vf); % MMA parameter, tunings are allowed 
%% ---------------------------Iteration-----------------------------
while (max(abs(obj./obj_old-1))>1e-3||abs(con)/vf>1e-2||Iter<=50)&&Iter<300
    %% ----------------------FE analysis ---------------------------
    FE_Phi = conv2(Phi,ones(2)/4,'valid'); % LS values on FE nodes
    Gauss_Phi = [(1-gx(:)).*(1-gy(:)),(1-gx(:)).*(1+gy(:)),(1+gx(:))...
        .*(1+gy(:)),(1+gx(:)).*(1-gy(:))]/4*FE_Phi(FEnod_loc);
    ke = Gauss_k*(1e-3*(Gauss_Phi<0)+(Gauss_Phi>=0)*1);
    K = sparse(Ipp(half),Jpp(half),ke(half));    K = K+K'-diag(diag(K));
    u = zeros(size(f));  u(freedof) = K(freedof,freedof)\f(freedof);
    obj_old = [obj_old(2:end);obj];   obj = f'*u; % Objective function
    con = sum(gweight'*(Gauss_Phi>=0))/4/(nx*ny)-vf; % Material volume
    figure(1); contourf(FE_Phi,[0 0]); colormap([0.6 0 0]); axis equal off;
    fprintf('It.: %4i  Vol_cons.: %6.3f Obj.: %6.3f\n',Iter,con,obj);
    %% -----------------------Sensitivity---------------------------
    dirac = 3/4/0.5*(1-Gauss_Phi.^2/0.5^2).*(abs(Gauss_Phi-0.5/2)<0.5/2);
    dcon(2:ny+1,2:nx+1) = reshape(gweight'*dirac/4/nx/ny,ny,nx);
    ku_e = reshape(Gauss_k*dirac*(1-1e-3).*u(repmat(FE_dof,8,1)),8,8*nx*ny);
    uku_e = reshape(sum(ku_e),8,nx*ny).*u(FE_dof);
    dobj(2:ny+1,2:nx+1) = -reshape(sum(uku_e),ny,nx);
    %% ---------------Optimization (MMA) & Evolution----------------
    [Vn0,~,~,~,~,~,~,~,~,low,upp] = mmasub(1,prod(LSgrid),Iter,Vn(:), ...
        -Vn_max,Vn_max,Vn_old1(:),Vn_old2(:),obj,dobj(:),con,dcon(:)', ...
        low,upp,1,0,cc(obj,con),0); % MMA, 09.2007 (a small change 08.2008)
    Vn_old2 = Vn_old1;  Vn_old1 = Vn; Vn = reshape(roundn(Vn0,-6),LSgrid);
    [Phi] = Evol_rein(Phi,Vn,Vn_max/1e3,120,Vn_max/10);  Iter = Iter+1;
end
%% -------------------------Subfunctions----------------------------
function [Phi] = Evol_rein(Phi,Vn,kesi,re_num,dt_re) % H-J & reinit.
[Grad_P,Grad_M,Grad2,~] = Phi_Diff(Phi);
Phi = Phi+(max(Vn,0).*Grad_P+min(Vn,0).*Grad_M)+kesi*Grad2;
[~ ,~,~,S] = Phi_Diff(Phi);
for i = 1:re_num % Reinitialization
    [Grad_P ,Grad_M,~,~] = Phi_Diff(Phi);
    Phi = Phi-dt_re/max(abs(S(:))).*((max(S,0).*Grad_P+min(S,0).*Grad_M)-S);
end
Phi = roundn(Phi,-6);
function [Grad_P,Grad_M,Grad2,S] = Phi_Diff(Phi) % LS Gradient (upwind)
DxL = (Phi-[ Phi(:,1)*2-Phi(:,2),Phi(:,1:end-1)]);
DxR = ([Phi(:,2:end),Phi(:,end)*2-Phi(:,end-1)]-Phi);
DyL = (Phi-[Phi(1,:)*2-Phi(2,:);Phi(1:end-1,:)]);
DyR = ([Phi(2:end,:);Phi(end,:)*2-Phi(end-1, :)] - Phi);
Grad_P = sqrt(max(DxL,0).^2+min(DxR,0).^2+max(DyL,0).^2+min(DyR,0).^2);
Grad_M = sqrt(min(DxL,0).^2+max(DxR,0).^2+min(DyL,0).^2+max(DyR,0).^2);
[Phixx,Phiyy] = deal(DxR-DxL,DyR-DyL);   Grad2 = Phixx+Phiyy;
S = Phi./(sqrt(Phi.^2+(((DxL+DxR)/2).^2+((DyL+DyR)/2).^2))+eps); % sign(Phi)