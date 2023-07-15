%%% A 51 line code for topology optimization with same size elements. Vittorio Latorre Dec. 2018 %%%
function [c,xnew,loop]=tosse(nelx,nely,volfrac,mu)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
U = zeros(2*(nely+1)*(nelx+1),1);
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% INITIALIZE ITERATION
vgam = 1;
loop = 0;
change = 1;
xPhys = ones(nely,nelx)*vgam;
%% START ITERATION
while change >   1e-16 &&loop<200
    loop = loop + 1;
    vgam = max(volfrac,vgam*mu);
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK);
    K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = (xPhys(:)+Emin)*E0.*(sum((U(edofMat)*KE).*U(edofMat),2));
    c = sum((Emin+ce*(1-Emin/E0)));
    %% SOLVE THE KNAPSACK PROBLEM
    [~,I]=sort(ce,1,'descend');
    xnew=zeros(nelx*nely,1);
    xnew(I(1:floor(vgam*nelx*nely)))=1;
    change = max(abs(xnew(:)-xPhys(:)));
    xPhys = reshape(xnew,[nely,nelx]);
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
end
%% PLOT DENSITIES
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;