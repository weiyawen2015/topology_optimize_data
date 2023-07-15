function VFLSM_3D(nx,ny,nz,vf,xnum) % A Matlab (2010a or later) code for VFLS (3D)
%% -------------------FE mesh & Level Set (LS) grid-----------------
nod_n = [ny+1,nx+1,nz+1]; FEnod_glob = reshape(1:prod(nod_n),nod_n); % Node
FEnod_loc = reshape(FEnod_glob(1:ny,1:nx,1:nz),1,nx*ny*nz); % Local node index
FEnod_loc = repmat(FEnod_loc,4,1)+repmat([0;1;ny+2;ny+1],1,nx*ny*nz);
FEnod_loc = [FEnod_loc;FEnod_loc+(nx+1)*(ny+1)]; % Node index
LSgrid = nod_n+1; [LSy,LSx,LSz] = ndgrid(-0.5:ny+0.5,-0.5:nx+0.5,-0.5:nz+0.5);
[fx,fy,fz,fixed] = deal(nx/2+1,ny/2+1,nz+1,[1:3,3*ny+3,3*(ny+1)*nx+[2,3,3*ny+3]]);
f = sparse(3*((ny+1)*(nx+1)*(fz-1)+(ny+1)*(fx-1)+fy),1,-1,3*prod(nod_n),1);
freedof = 1:3*prod(nod_n); freedof(fixed) = []; % Free node
fixed_corn = [1,ny+2,(ny+2)*(nx+1)+1,(ny+2)*(nx+2)];
FE_dof = 3*FEnod_loc(repmat(1:8,3,1),:)-repmat([2;1;0],8,nx*ny*nz); % DOF index
[Ipp,Jpp] = deal(FE_dof(:,repmat(1:nx*ny*nz,24,1)),FE_dof(repmat(1:24,24,1),:));
half = (Ipp(:)<=Jpp(:)); % Preparation for assembling stiffness matrix
%% --------------------------Gauss points---------------------------
gs = [(2*70^0.5)/63+5/9,5/9-(2*70^0.5)/63].^0.5;   gauss = [-gs,0,gs];
gw = 2./(1-gauss.^2)./((315*gauss.^4)/8-(105*gauss.^2)/4+15/8).^2;
[gy,gx,gz] = ndgrid(gauss,gauss,gauss); [gwy,gwx,gwz] = ndgrid(gw,gw,gw);
gweight = gwx(:).*gwy(:).*gwz(:); gnum = length(gweight); % Weighting & number
NGauss = [(1-gx(:)).*(1-gy(:)),(1-gx(:)).*(1+gy(:)),...
    (1+gx(:)).*(1+gy(:)),(1+gx(:)).*(1-gy(:))]/8;
NGauss = [NGauss.*(1-repmat(gz(:),1,4)),NGauss.*(1+repmat(gz(:),1,4))];
Dmat = (1-0.3)/(1+0.3)/(1-2*0.3)*diag(kron([1,(1-0.3*2)/2*(1-0.3)],[1,1,1]));
Dmat([2,3,7,9,13,14]) = 0.3/(1-0.3);
Bgsx = @(i)kron([1-gz(i),1+gz(i)],[-1+gy(i),-gy(i)-1,1+gy(i),1-gy(i)]/4);
Bgsy = @(i)kron([1-gz(i),1+gz(i)],[-1+gx(i),1-gx(i),gx(i)+1,-gx(i)-1]/4);
Bgsz = @(i)kron([-1,1],[(1-gx(i))*(1-gy(i)),(1-gx(i))*(1+gy(i)),...
    (1+gx(i))*(1+gy(i)),(1+gx(i)).*(1-gy(i))]/4);
Bgs = @(i)[reshape([Bgsx(i);zeros(2,8)],1,24);reshape([zeros(1,8);Bgsy(i);...
    zeros(1,8)],1,24);reshape([zeros(2,8);Bgsz(i)],1,24);reshape([Bgsy(i);...
    Bgsx(i);zeros(1,8)],1,24);reshape([zeros(1,8);Bgsz(i);Bgsy(i)],1,24);...
    reshape([Bgsz(i);zeros(1,8);Bgsx(i)],1,24)];
BDB = @(i)reshape(Bgs(i)'*Dmat*Bgs(i),24^2,1)*gweight(i)/8;
Gauss_k = cell2mat(arrayfun(@(i)BDB(i),1:gnum,'UniformOutput',false));
%% -------------------------Initialization--------------------------
[ynum,znum,cr] = deal(floor(xnum/nx*ny),floor(xnum/nx*nz),nx/xnum/4);
[cdx,cdy,cdz] = deal(nx/xnum/2,ny/ynum/2,nz/znum/2);
[cy,cx,cz] = ndgrid((0:2:ynum*2+2)*cdy,(1:2:xnum*2+3)*cdx,(0:2:znum*2+2)*cdz);
[cx,cy,cz] = deal([cx;cx-cdx],[cy;cy+cdy],[cz;cz]);
[cx,cy,cz] = deal([cx(:);cx(:)-cdx],[cy(:);cy(:)],[cz(:);cz(:)+cdz]);
cNAN = max(abs([cx-fx,cy-fy,cz-fz])-1,[],2)<=cr;
[cx(cNAN),cy(cNAN),cz(cNAN)] = deal(inf);
[cx,xt] = ndgrid(cx,LSx);  [cy,yt] = ndgrid(cy,LSy); [cz,zt] = ndgrid(cz,LSz);
Phi = sqrt((xt-cx).^2+(yt-cy).^2+(zt-cz).^2)-cr; % Initial LS function
Phi = min([Phi;[LSx(:),nx-LSx(:),LSy(:),ny-LSy(:),LSy(:),nz-LSz(:),LSz(:)]'+0.5]);
Phi(fixed_corn) = 0.5;    Phi = reshape(Phi,LSgrid); % LS values on FE nodes
Iter = 0;  obj_old = inf(3,1);  obj = 0;  con = 1;
[Vn,Vn_old1,Vn_old2,dobj,dcon] = deal(zeros(LSgrid)); % Des.var.& sensi.
[Vn_max,upp,low] = deal(0.1,0.1,-0.1); % Des. var. bounds & asymptote in MMA 
cc = @(obj,con)0.75*obj/max((con+vf),vf); % MMA parameter, tunnings are allowed
%% ---------------------------Iteration-----------------------------
while (max(abs(obj./obj_old-1))>1e-3||abs(con)/vf>2e-2||Iter<=50)&&Iter<300
    %% ----------------------FE analysis ---------------------------
    FE_Phi = convn(Phi,ones(2,2,2)/8,'valid'); % LS values on FE nodes
    Gauss_Phi = NGauss*FE_Phi(FEnod_loc);
    ke = Gauss_k*(1e-3*(Gauss_Phi<0)+(Gauss_Phi>=0)*1);
    K = sparse(Ipp(half),Jpp(half),ke(half));    K = K+K'-diag(diag(K));
    u = zeros(size(f));  u(freedof) = K(freedof,freedof)\f(freedof);
    obj_old = [obj_old(2:end);obj];  obj = f'*u; % Objective function
    con = sum(gweight'*(Gauss_Phi>=0))/8/(nx*ny*nz)-vf; % Material volume
    cla(); patch(isosurface(FE_Phi,0),'FaceColor','r','EdgeColor','none');
    patch(isocaps(FE_Phi,0),'FaceColor','r','EdgeColor','none');
    camlight right; view ([145,25]); axis equal off; lighting gouraud; drawnow;
    fprintf('It.: %4i  Vol_cons.: %7.4f Obj.: %7.4f\n',Iter,con,obj);
    %% -----------------------Sensitivity---------------------------
    dirac = 3/4/0.5*(1-Gauss_Phi.^2/0.5^2).*(abs(Gauss_Phi-0.5/2)<0.5/2);
    dcon(2:ny+1,2:nx+1,2:nz+1) = reshape(gweight'*dirac/8/nx/ny/nz,ny,nx,nz);
    ku_e = reshape(Gauss_k*dirac*(1-1e-3).*u(repmat(FE_dof,24,1)),24,24*nx*ny*nz);
    uku_e = reshape(sum(ku_e),24,nx*ny*nz).*u(FE_dof);
    dobj(2:ny+1,2:nx+1,2:nz+1) = -reshape(sum(uku_e),ny,nx,nz);
    %% ---------------Optimization (MMA) & Evolution----------------
    [Vn0,~,~,~,~,~,~,~,~,low,upp] = mmasub(1,prod(LSgrid),Iter,Vn(:),...
        -Vn_max,Vn_max,Vn_old1(:),Vn_old2(:),obj,dobj(:),con,dcon(:)',...
        low,upp,1,0,cc(obj,con),0); % MMA, 09.2007 (a small change 08.2008)
    Vn_old2 = Vn_old1; Vn_old1 = Vn; Vn = reshape(roundn(Vn0,-6),LSgrid);
    [Phi] = Evol_rein(Phi,Vn,Vn_max/1e3,120,Vn_max/10);  Iter = Iter+1;
end
%% -------------------------Subfunctions----------------------------
function [Phi] = Evol_rein(Phi,Vn,kesi,re_num,dt_re) % H-J & reinit.
[Grad_P,Grad_M,Grad2,~] = Phi_Diff(Phi);
Phi = Phi+((max(Vn,0).*Grad_P+min(Vn,0).*Grad_M)+kesi*Grad2);
[~ ,~,~,S] = Phi_Diff(Phi); % sign(Phi)
for i = 1 : re_num % Reinitialization
    [Grad_P,Grad_M,~,~] = Phi_Diff(Phi);
    Phi = Phi-dt_re/max(abs(S(:))).*((max(S,0).*Grad_P+min(S,0).*Grad_M)-S);
end
Phi = roundn(Phi,-6);
function [Grad_P,Grad_M,Grad2,S] = Phi_Diff(Phi) % LS Gradient (upwind)
DxL = (Phi-cat(2,Phi(:,1,:)*2-Phi(:,2,:),Phi(:,1:end-1,:)));
DxR = (cat(2,Phi(:,2:end,:),Phi(:,end,:)*2-Phi(:,end-1,:))-Phi);
DyL = (Phi-cat(1,Phi(1,:,:)*2-Phi(2,:,:),Phi(1:end-1,:,:)));
DyR = (cat(1,Phi(2:end,:,:),Phi(end,:,:)*2-Phi(end-1, :,:)) - Phi);
DzL = (Phi-cat(3,Phi(:,:,1)*2-Phi(:,:,2),Phi(:,:,1:end-1)));
DzR = (cat(3,Phi(:,:,2:end),Phi(:,:,end)*2-Phi(:,:,end-1)) - Phi);
Grad_P = sqrt(max(DxL,0).^2+min(DxR,0).^2+max(DyL,0).^2+min(DyR,0).^2 ...
    +max(DzL,0).^2+min(DzR,0).^2);
Grad_M = sqrt(min(DxL,0).^2+max(DxR,0).^2+min(DyL,0).^2+max(DyR,0).^2 ...
    +min(DzL,0).^2+max(DzR,0).^2);
[Phixx,Phiyy,Phizz] = deal(DxR-DxL,DyR-DyL,DzR-DzL); Grad2 = Phixx+Phiyy+Phizz;
S = Phi./(sqrt(Phi.^2+(((DxL+DxR)/2).^2+((DyL+DyR)/2).^2+((DzL+DzR)/2).^2))+eps);