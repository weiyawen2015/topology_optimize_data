%%%% CODE of Combining BESO with RL BY H.B.SUN and L.MA in 2019 (Geometrical nonlinerity) %%%%
clc
clear all
close all
%========================== preprocess
dbstop if error

nelx = 60; nely = 40;
lx = 60; ly = 40;
volfrac = 0.5;
E = 1e3;
nu = 0.3;
penal = 3;
er = 0.03;
[nodes,eles,eleNum,nodeNum] = GenerateMesh(lx,ly,nelx,nely);

rmin = 3;
Pfinal = 600/10; % the given pressure
totalstep = 10;
% deltaP = Pfinal/totalstep;

% ========================  Preparation FE analysis ======================
fixeddofs=[1:2*(nely+1)]; %left
loaddofs = 2*(nelx+1)*(nely+1)-nely;
alldofs = 1:2*(nely+1)*(nelx+1);
alldofNum = 2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
pic=zeros(100,1);R=zeros(100,nelx*nely);V=zeros(100,nelx*nely);gamma=0.9;winx=0;winy=0;
Nsubopt=20;unfixvol=1-volfrac/2; emax=0.1;w=0.005;
WC=zeros(Nsubopt+1,1);ITER=zeros(Nsubopt+1,1);WC_difference=zeros(Nsubopt,1);IOU=zeros(Nsubopt,1);

F = sparse(alldofNum,1);
uhat = zeros(alldofNum,1);  % for sensitivity: K_T*u = F
uout = zeros(alldofNum,1);  % for real displacement
deltadisp = zeros(alldofNum,1);

% save dof of every element
edofMat = zeros(eleNum,8);
for i = 1:eleNum
    node1 = eles(i,1); 
    node2 = eles(i,2); 
    node3 = eles(i,3); 
    node4 = eles(i,4);
    edofMat(i,:) = [2*node1-1,2*node1, 2*node2-1,2*node2,2*node3-1,2*node3,2*node4-1,2*node4];
end
% ======================================================================

[H,Hs] = Filter(nelx,nely,rmin);

dc=zeros(nely,nelx);
dv=zeros(nely,nelx);

% =========================Parameter of MMA =============================
xval  = ones(1,nely,nelx);
m     = 1;                	 % The number of general constraints.
xmin  = zeros(nely*nelx,1);
xmax  = ones(nely*nelx,1);
xold1 = xval(:);             % xval, one iteration ago (provided that iter>1).
xold2 = xval(:);             % xval, two iterations ago (provided that iter>2).
low   = xmin;
upp   = xmax;
a0    = 1;                	 % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       	 % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 1000*ones(m,1);   	 % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);
% =======================================================================

% pre caculate to speed up
[wt,shape,dshapedxi,dshapedeta] = shapeFunction();
[detjacob,invjacob] = Jacobian(nodes,eles,dshapedxi,dshapedeta,eleNum);
[dshapedx,dshapedy,G] = dNdx_dNdy_G(invjacob,dshapedxi,dshapedeta,eleNum);
B0 = kinematicStiffnessLinear(dshapedx,dshapedy,eleNum);

loop=0; change=1; maxiter=100;
f0val = zeros(maxiter,1); Wc = zeros(maxiter,1);
vol = 1; 
while (change > 0.01)&&(loop<maxiter)
    loop=loop+1;
    vol = max(vol*(1-er),volfrac);
    if loop >1; 
        olddc = dc; 
    end
	% Newton raphson iteration to update displacement uout
    P = 0; u = 0;
    
    uout(:) = 0;
    uhat(:) = 0;
    deltadisp(:) = 0;
    for loadstep = 1:totalstep+1
        % For first iteration the increment in load is very small, as it is imp to get first solution correct.
        if loadstep==1
            deltaP = 0.01*Pfinal;
        else  
            deltaP = Pfinal/totalstep;
        end
        P = P+deltaP;
        F(loaddofs,1) = P;
                
        % until residual approximately zero
        for updateintforce = 1:100
            [stiffness,tangentstiffness] = FEA(eleNum,uout,alldofNum,E,nu,edofMat,...
                fixeddofs,xval(:),penal,wt,dshapedx,dshapedy,detjacob,G,B0);
            
            residual = stiffness(freedofs,freedofs)*uout(freedofs) - F(freedofs,1);
            
            deltadisp(freedofs) = -tangentstiffness(freedofs,freedofs)\residual;
            uout = uout+deltadisp;
            
            residualnorm = norm(residual);
            if residualnorm<10^-3
                break; % jump out of "for"
            end
            
        end   
        Wc(loop) = Wc(loop) + 0.5*(uout(loaddofs,1)+u)*deltaP;
        u = uout(loaddofs,1);
    end
    uhat(freedofs) = tangentstiffness(freedofs,freedofs) \ F(freedofs);
    [c,dc] = SS(eleNum,uout,E,nu,...
        xval(:),penal,wt,dshapedx,dshapedy,detjacob,edofMat,uhat,B0,G);

    dv = ones(nely,nelx)/(volfrac*eleNum);
	
	% filter to remove checkboard
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);

    % BESO optimization to update density xval 
    % STABLIZATION OF EVOLUTIONARY PROCESS
    if loop > 1; 
        dc = (dc+olddc)/2.; 
    end
    
    [xval]=ADDDEL(nelx,nely,vol,dc);
    if (vol-volfrac)>0.01
        qi=loop+5;
    end    
    if qi>=loop
        R(loop,:)=reshape(xval,1,nelx*nely);
    end
    volnext(loop)=mean(xval(:));
    
    f0val(loop) = sum(c);
    if loop>10;
        change=abs(sum(f0val(loop-9:loop-5))-sum(f0val(loop-4:loop)))/sum(f0val(loop-4:loop));
    end
    
	% output
    fprintf(' It.:%4i Obj.:%6.3f W.:%6.3f Vol.:%6.4f ch.:%6.4f \n',loop,f0val(loop),Wc(loop),vol,change);
    colormap(gray); imagesc(1-xval); caxis([0 1]); axis equal; axis off; drawnow;
   
end
W0=Wc(loop);
WC(1)=W0;ITER(1)=loop;
xsubopt=reshape(xval,1,nelx*nely);nsubopt=[];N=0;

%% RL section
flag = 0;
for z=1:Nsubopt
if Wc(loop)/W0-1<0.2 & flag == 0
    N=N+1;
    for j=qi-1:-1:1
        R(j,:)=R(j,:)+gamma*R(j+1,:);
    end
    V(1:qi,:)=V(1:qi,:)+1/N*(R(1:qi,:)*(2-Wc(loop)/W0)-V(1:qi,:));
end
flag = 0;
% INITIALIZE
xval(1:nely,1:nelx) = 1.;
vol=1; loop = 0; change = 1.; 
R=zeros(100,nelx*nely);
f0val = zeros(maxiter,1); Wc = zeros(maxiter,1);
% START iTH ITERATION
while (change > 0.01)&&(loop<maxiter)||(vol>volfrac)
    loop=loop+1;
    vol = max(vol*(1-er),volfrac);
    if loop >1; 
        olddc = dc; 
    end
	% Newton raphson iteration to update displacement uout
    P = 0;  u = 0;
    uout(:) = 0;
    uhat(:) = 0;
    deltadisp(:) = 0;
    for loadstep = 1:totalstep+1
        % For first iteration the increment in load is very small, as it is imp to get first solution correct.
        if loadstep==1
            deltaP = 0.01;
        else  
            deltaP = Pfinal/totalstep;
        end
        P = P+deltaP;
        F(loaddofs,1) = P;
                
        % until residual approximately zero
        for updateintforce = 1:100
            [stiffness,tangentstiffness] = FEA(eleNum,uout,alldofNum,E,nu,edofMat,...
                fixeddofs,xval(:),penal,wt,dshapedx,dshapedy,detjacob,G,B0);
            
            residual = stiffness(freedofs,freedofs)*uout(freedofs) - F(freedofs,1);
            deltadisp(freedofs) = -tangentstiffness(freedofs,freedofs)\residual;
            uout = uout+deltadisp;
            
            residualnorm = norm(residual);
            if residualnorm<10^-3
                break; % jump out of "for"
            end
            
        end 
        if flag == 1
            break;
        end
        Wc(loop) = Wc(loop) + 0.5*(uout(loaddofs,1)+u)*deltaP;
        u = uout(loaddofs,1);
    end
    if flag == 1
        break;
    end
    uhat(freedofs) = tangentstiffness(freedofs,freedofs) \ F(freedofs);
    [c,dc] = SS(eleNum,uout,E,nu,...
        xval(:),penal,wt,dshapedx,dshapedy,detjacob,edofMat,uhat,B0,G);

    dv = ones(nely,nelx)/(volfrac*eleNum);
	
	% filter to remove checkboard
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);

    % BESO optimization to update density xval 
    % STABLIZATION OF EVOLUTIONARY PROCESS
    if loop > 1; 
        dc = (dc+olddc)/2.; 
    end
    
    e=emax+(0-emax)*min(1,1/(qi-5)*loop)^3;
    if (vol-volfrac)>0.01;
        [xval] = RL_ADDDEL(loop,qi,nelx,nely,unfixvol,volnext(loop),V(loop,:),dc,e,winx,winy,w);
    else
        [xval] = ADDDEL(nelx,nely,vol,dc);
    end
    if qi>=loop
        R(loop,:)=reshape(xval,1,nelx*nely);
    end
	% output
    f0val(loop) = sum(c);
    if loop>10;
        change=abs(sum(f0val(loop-9:loop-5))-sum(f0val(loop-4:loop)))/sum(f0val(loop-4:loop));
    end
    fprintf(' It.:%4i Obj.:%6.3f W.:%6.3f Vol.:%6.4f ch.:%6.4f \n',loop,f0val(loop),Wc(loop),vol,change);
    figure(z+1);
    colormap(gray); imagesc(1-xval); caxis([0 1]); axis equal; axis off; drawnow;

end
if flag == 1
    continue;
end
WC(z+1)=Wc(loop);ITER(z+1)=loop;
WC_difference(z)=WC(z+1)/WC(1)-1;
if WC_difference(z)<0.2
    xval=reshape(xval,1,nelx*nely);
    for i=1:size(xsubopt,1)
        IOU(z)=max(IOU(z),length(find((xval+xsubopt(i,:))>1.5))/length(find((xval+xsubopt(i,:))>0.5)));
        if IOU(z)>0.9
            break;
        end
        if i==size(xsubopt,1)
            xsubopt=[xsubopt;xval];
            nsubopt=[nsubopt;z+1];
        end
    end
end
if length(nsubopt)==10
    break
end
end
% for z=1:length(nsubopt)
%     figure(1+z);
%     xval=reshape(xsubopt(1+z,:),nely,nelx);
%     colormap(gray); imagesc(-xval); axis equal; axis tight; axis off;
% end

%%%%%%%%%% BESO UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=ADDDEL(nelx,nely,vol,dc)
dc=reshape(dc,nely,nelx);
l1 = min(min(dc)); l2 = max(max(dc));
while ((l2-l1)/l2 > 1.0e-5)
    th = (l1+l2)/2.0;
    x = max(0.001,sign(dc-th));
    if sum(sum(x))-vol*(nelx*nely) > 0;
        l1 = th;
    else
        l2 = th;
    end
end
end

%%%%%%%%%% RL UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=RL_ADDDEL(i,qi,nelx,nely,unfixvol,volnext,V,dc,e,winx,winy,w)
% % symmetric constraint
% if sym==1
%     V=reshape(V,nely,nelx);
%     nely=nely/2;
%     V=V(1:nely,:);
%     V=reshape(V,1,nelx*nely);
%     dc=dc(1:nely,:);
% end 
% determine V by both V and dc
[dc0,index]=sort(dc);
Vnow=min(V)+(dc-dc0(1))/(dc0(unfixvol*nelx*nely)-dc0(1))*(max(V)-min(V));
V=w*V+(1-w)*Vnow;

% delete 1-e elements
[V0,index]=sort(V);
base=floor((1-e)*(1-volnext)*(nelx*nely));
x=ones(1,nely*nelx);
V_critical=V0(base); % RANDOM CHOICE
index_critical=find(V0==V_critical);
base2=index_critical(1);
index(index_critical)= index(index_critical(1)-1+randperm(length(index_critical)));
x(index(1:base))=0.001;
x=reshape(x,nely,nelx);
% delete e elements by window
m=floor(e*(1-volnext)*(nelx*nely)/(2*winx+1)/(2*winy+1));
if m>0
n=unfixvol*nelx*nely-base;
p=randperm(n);
for i=1:m
    ind=index(base+p(i));
    indx=ceil(ind/nely);
    indy=mod(ind,nely);
    if indy==0
        indy=nely;
    end
    x(indy,max(1,indx-winx):min(indx+winx,nelx))=0.001;
    if winy>0
        for j=1:winy
            x(max(1,indy-j),max(1,indx-winx):min(indx+winx,nelx))=0.001;
            x(min(nely,indy+j),max(1,indx-winx):min(indx+winx,nelx))=0.001;
        end
    end
end
end
% delete rest elements
n=floor((mean(x(:))-volnext)*(nelx*nely));
if n>0
    ind2=find(x(index(1:unfixvol*nelx*nely))==1);
    p=randperm(length(ind2));
    x(index(ind2(p(1:n))))=0.001;
end
% if sym==1
%     x=[x;flipud(x)];
% end 
end
