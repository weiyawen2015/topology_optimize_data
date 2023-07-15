function [stiffness,tangentstiffness] = FEA(eleNum,displacement,alldofNum,E,nu,edofMat,...
    fixeddofs,x,penal,wt,dshapedx,dshapedy,detjacob,G,B0)
%--------------------------------------------------------------------------
%  Computation of element matrices and vectors and their assembly
%--------------------------------------------------------------------------

tangentstiffness = zeros(alldofNum,alldofNum) ;      % System tangent stiffness matrix
stiffness = zeros(alldofNum,alldofNum);              % system stiffness matrix


D = E/(1-nu^2)*[1  nu 0; nu  1  0; 0  0  (1-nu)/2];  % material matrix

for i = 1:eleNum          
    edof = edofMat(i,:);
    u = displacement(edof(1:2:end));
    v = displacement(edof(2:2:end)); 
    
    kb = zeros(8,8);              % initialization of element nonlinear stiffness matrix
    kl = zeros(8,8);              % initialization of element nonlinear stiffness matrix
    kll = zeros(8,8);             % initialization of element nonlinear stiffness matrix
    ksigma=zeros(8,8);            % initialization of element initial stress matrix
    
    % Gauss integration
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
        kl = kl+(0.5*B0i'*D*BLi+BLi'*D*B0i+0.5*BLi'*D*BLi)*wt(j)*detjacob(i,j);
        kll= kll+(B0i'*D*BLi+BLi'*D*B0i+BLi'*D*BLi)*wt(j)*detjacob(i,j);
    end
        
    ke = kb+kl;
    kt = kb+kll+ksigma;
    
    % assemble
    stiffness(edof,edof) = stiffness(edof,edof)+ke*x(i)^penal;
    tangentstiffness(edof,edof) = tangentstiffness(edof,edof)+kt*x(i)^penal;
end

% constrains
% for i = 1:length(fixeddofs)
%     for j = 1:alldofNum
%         stiffness(fixeddofs(i),j) = 0;
%         stiffness(j,fixeddofs(i)) = 0;
%         tangentstiffness(fixeddofs(i),j) = 0;
%         tangentstiffness(j,fixeddofs(i)) = 0;
%     end
%     stiffness(fixeddofs(i),fixeddofs(i)) = 1;
%     tangentstiffness(fixeddofs(i),fixeddofs(i)) = 1;
% end