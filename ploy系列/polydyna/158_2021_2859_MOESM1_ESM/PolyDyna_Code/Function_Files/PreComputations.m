%--------------------------- PreComputations -----------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    %
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
function fem = PreComputations(fem)
%Element areas
fem.ElemArea = zeros(fem.NElem,1); 
for el = 1:fem.NElem 
  vx = fem.Node(fem.Element{el},1); vy = fem.Node(fem.Element{el},2);
  temp = vx.*vy([2:end 1])-vy.*vx([2:end 1]);
  fem.ElemArea(el,1) = 0.5*sum(temp);
end
%Shape functions for FE analysis
fem = TabShapeFnc(fem);
%Compute triplets used to assemble stiffness matrix and load vector
fem.ElemNDof = 2*cellfun(@length,fem.Element); % # of DOFs per element
fem.i = zeros(sum(fem.ElemNDof.^2),1); fem.j=fem.i; fem.e=fem.i; 
fem.k0 = fem.i; fem.m0 = fem.i; fem.Fa0 = zeros(sum(fem.ElemNDof),1); 
fem.eDof = fem.Fa0; %fem.DofE = fem.eDof; fem.Fa0=fem.eDof;
fem.iK0=fem.i; fem.jK0=fem.iK0;
indexK = 0; indexF = 0; IniI=0; IniJ=0;
for el = 1:fem.NElem  
  NDof = fem.ElemNDof(el); NKE = NDof^2;
  eDof = reshape([2*fem.Element{el}-1;2*fem.Element{el}],NDof,1); 
  %Triplets for stiffness matrix
  I = repmat(eDof ,1,NDof); J = I';
  fem.i(indexK+1:indexK+NKE) = I(:);
  fem.j(indexK+1:indexK+NKE) = J(:); 
  fem.e(indexK+1:indexK+NKE) = el;
  fem.DofE(indexF+1:indexF+NDof) = el; 
  fem.eDof(indexF+1:indexF+NDof) = eDof;
  %Local stiffness and mass matrices
  if ~fem.Reg || el==1,  [Ke, Me] = LocalKM(fem,fem.Element{el}); end
  fem.k0(indexK+1:indexK+NKE) = Ke(:);
  fem.m0(indexK+1:indexK+NKE) = Me(:);
  fem.Fa0(indexF+1:indexF+NDof) = Me*ones(NDof,1);
  %Element DOFs 
  I = repmat((1:NDof)',1,NDof); J=I';
  fem.iK0(indexK+1:indexK+NKE) = I(:)+IniI;
  fem.jK0(indexK+1:indexK+NKE) = J(:)+IniJ;
  IniI = IniI+NDof;
  IniJ = IniJ+NDof;
  indexK = indexK+NDof^2; indexF = indexF+NDof; 
end
fem.k0 = fem.Thickness.*fem.k0;
fem.m0 = fem.Thickness.*fem.rho.*fem.m0;
if ~isempty(fem.ag), fem.Fa0 = -fem.Thickness.*fem.Fa0.*fem.ag;
else, fem.Fa0 = fem.Fa0.*zeros(1,fem.NStep+1); end
%Compute external load vector and initialize inertia force vector
NLoad = size(fem.Load,1);
fem.Fext = zeros(2*fem.NNode,fem.NStep+1);  %External load vector
fem.Fa = 0.*fem.Fext;                       %Inertia force vector
if NLoad>0
 fem.Fext(2*fem.Load(:,1)-1,:) = fem.Load(:,2:2:end);  %x-crdnt
 fem.Fext(2*fem.Load(:,1),:)   = fem.Load(:,3:2:end);  %y-crdnt
end
%Obtain fixed DOFs and free DOFs
NSupp = size(fem.Supp,1);
FixedDofs = [fem.Supp(1:NSupp,2).*(2*fem.Supp(1:NSupp,1)-1);
             fem.Supp(1:NSupp,3).*(2*fem.Supp(1:NSupp,1))];
FixedDofs = FixedDofs(FixedDofs>0);
AllDofs   = 1:2*fem.NNode;
fem.FreeDofs = setdiff(AllDofs,FixedDofs);
%-------------------------------------- ELEMENT STIFFNESS AND MASS MATRICES
function [Ke,Me] = LocalKM(fem,eNode)
D = fem.E0/(1-fem.Nu0^2)*[1 fem.Nu0 0;fem.Nu0 1 0;0 0 (1-fem.Nu0)/2];
nn = length(eNode); Ke = zeros(2*nn,2*nn); Me = Ke;
W = fem.ShapeFnc{nn}.W;
B = zeros(3,2*nn); N = zeros(2,2*nn); 
for q = 1:length(W)  %Quadrature loop
  N(1,1:2:end) = fem.ShapeFnc{nn}.N(:,:,q); N(2,2:2:end) = N(1,1:2:end);
  dNdxi = fem.ShapeFnc{nn}.dNdxi(:,:,q);
  J0 = fem.Node(eNode,:)'*dNdxi; WdJ0 = W(q)*det(J0);
  dNdX = dNdxi/J0;
  B(1,1:2:2*nn) = dNdX(:,1)'; B(2,2:2:2*nn) = dNdX(:,2)'; 
  B(3,1:2:2*nn) = dNdX(:,2)'; B(3,2:2:2*nn) = dNdX(:,1)';
  Ke = Ke+B'*D*B*WdJ0; %Stiffness matrix assuming unit thickness
  Me = Me+N'*N*WdJ0;   %Mass matrix assuming unit density
end
%------------------------------------------------- TABULATE SHAPE FUNCTIONS
function fem = TabShapeFnc(fem)
ElemNNode = cellfun(@length,fem.Element); %Number of nodes per element
fem.ShapeFnc = cell(max(ElemNNode),1);
for nn = min(ElemNNode):max(ElemNNode)
  [W,Q] = PolyQuad(nn);
  fem.ShapeFnc{nn}.W = W;
  fem.ShapeFnc{nn}.N = zeros(nn,1,size(W,1));
  fem.ShapeFnc{nn}.dNdxi = zeros(nn,2,size(W,1));
  for q = 1:size(W,1)
    [N,dNdxi] = PolyShapeFnc(nn,Q(q,:));
    fem.ShapeFnc{nn}.N(:,:,q) = N;
    fem.ShapeFnc{nn}.dNdxi(:,:,q) = dNdxi;
  end
end
%------------------------------------------------ POLYGONAL SHAPE FUNCTIONS
function [N,dNdxi] = PolyShapeFnc(nn,xi)
N=zeros(nn,1); alpha=zeros(nn,1); dNdxi=zeros(nn,2); dalpha=zeros(nn,2);
sum_alpha=0.0; sum_dalpha=zeros(1,2); A=zeros(nn,1); dA=zeros(nn,2);
[p,Tri] = PolyTrnglt(nn,xi);
for i=1:nn
  sctr = Tri(i,:); pT = p(sctr,:);
  A(i) = 1/2*det([pT,ones(3,1)]);
  dA(i,1) = 1/2*(pT(3,2)-pT(2,2));
  dA(i,2) = 1/2*(pT(2,1)-pT(3,1));
end
A=[A(nn,:);A]; dA=[dA(nn,:);dA];
for i=1:nn
  alpha(i) = 1/(A(i)*A(i+1));
  dalpha(i,1) = -alpha(i)*(dA(i,1)/A(i)+dA(i+1,1)/A(i+1));
  dalpha(i,2) = -alpha(i)*(dA(i,2)/A(i)+dA(i+1,2)/A(i+1));
  sum_alpha = sum_alpha + alpha(i);
  sum_dalpha(1:2) = sum_dalpha(1:2)+dalpha(i,1:2);
end
for i=1:nn
  N(i) = alpha(i)/sum_alpha;
  dNdxi(i,1:2) = (dalpha(i,1:2)-N(i)*sum_dalpha(1:2))/sum_alpha;
end
%---------------------------------------------------- POLYGON TRIANGULATION
function [p,Tri] = PolyTrnglt(nn,xi)
p = [cos(2*pi*((1:nn))/nn); sin(2*pi*((1:nn))/nn)]';
p = [p; xi];
Tri = zeros(nn,3); Tri(1:nn,1)=nn+1;
Tri(1:nn,2)=1:nn; Tri(1:nn,3)=2:nn+1; Tri(nn,3)=1;
%----------------------------------------------------- POLYGONAL QUADRATURE
function [weight,point] = PolyQuad(nn)
[W,Q]= TriQuad;                  %Integration pnts & wgts for ref. triangle
[p,Tri] = PolyTrnglt(nn,[0 0]);  %Triangulate from origin
point=zeros(nn*length(W),2); weight=zeros(nn*length(W),1);
for k=1:nn
  sctr = Tri(k,:);
  for q=1:length(W)
    [N,dNds] = TriShapeFnc(Q(q,:));  %Compute shape functions
    J0 = p(sctr,:)'*dNds;
    l = (k-1)*length(W) + q;
    point(l,:) = N'*p(sctr,:);
    weight(l) = det(J0)*W(q);
  end                                 
end
%---------------------------------------------------- TRIANGULAR QUADRATURE
function [weight,point] = TriQuad
point=[1/6,1/6;2/3,1/6;1/6,2/3]; weight=[1/6,1/6,1/6];   
%----------------------------------------------- TRIANGULAR SHAPE FUNCTIONS
function [N,dNds] = TriShapeFnc(s)
N=[1-s(1)-s(2);s(1);s(2)]; dNds=[-1,-1;1,0;0,1];
%-------------------------------------------------------------------------%