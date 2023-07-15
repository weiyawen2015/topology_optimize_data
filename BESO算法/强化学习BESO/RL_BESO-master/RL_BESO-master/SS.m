function [c,dc] = SS(eleNum,displacement,E,nu,...
    x,penal,wt,dshapedx,dshapedy,detjacob,edofMat,uhat,B0,G)
%--------------------------------------------------------------------------
%  Computation of element compliance and derivatives and their assembly
%--------------------------------------------------------------------------

c = zeros(eleNum,1);  % element compliance
dc = zeros(eleNum,1); % element derivative of compliance

D = E/(1-nu^2)*[1  nu 0; nu  1  0; 0  0  (1-nu)/2];  % material matrix

for i = 1:eleNum
    edof = edofMat(i,:);
    u = displacement(edof(1:2:end));
    v = displacement(edof(2:2:end));
    
    kb = zeros(8,8);              
    kll = zeros(8,8);             
    ksigma=zeros(8,8);            
    
    for j = 1:4
        B0i = reshape(B0(i,(j-1)*3*8+1:j*3*8),3,8);
        Gi = reshape(G(i,(j-1)*4*8+1:j*4*8),4,8);
        
        dNdx = dshapedx(i,(j-1)*4+1:j*4);
        dNdy = dshapedy(i,(j-1)*4+1:j*4);
        dudx = dNdx*u;
        dudy = dNdy*u;
        dvdx = dNdx*v;
        dvdy = dNdy*v;
        
        straini = [dudx+dudx^2/2+dvdx^2/2;
            dvdy+dudy^2/2+dvdy^2/2;
            dudy+dvdx+dudx*dudy+dvdx*dvdy];
        
        c(i) = c(i) + straini'*D*straini*x(i)^penal*wt(j)*detjacob(i,j);
        
        stressi = D*straini;
        Si = [stressi(1)*eye(2),stressi(3)*eye(2);
            stressi(3)*eye(2),stressi(2)*eye(2)];
        
        BLi = zeros(3,8);
        BLi(1,1:2:end) = dudx*dNdx;
        BLi(1,2:2:end) = dvdx*dNdx;
        BLi(2,1:2:end) = dudy*dNdy;
        BLi(2,2:2:end) = dvdy*dNdy;
        BLi(3,1:2:end) = dudy*dNdx+dudx*dNdy;
        BLi(3,2:2:end) = dvdy*dNdx+dvdx*dNdy;
        
        ksigma = ksigma+Gi'*Si*Gi*wt(j)*detjacob(i,j);
        kb = kb+B0i'*D*B0i*wt(j)*detjacob(i,j);
        kll= kll+(B0i'*D*BLi+BLi'*D*B0i+BLi'*D*BLi)*wt(j)*detjacob(i,j);
    end
    kt = kb+kll+ksigma;
    dc(i) = uhat(edof)'*kt'*uhat(edof)*x(i)^(penal-1)*penal;
end