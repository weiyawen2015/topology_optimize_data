%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function fem = preComputations(fem)
%Element areas and von Mises matrix
fem.ElemArea = zeros(fem.NElem,1); 
for el = 1:fem.NElem 
  vx = fem.Node(fem.Element{el},1); vy = fem.Node(fem.Element{el},2);
  temp = vx.*vy([2:end 1])-vy.*vx([2:end 1]);
  fem.ElemArea(el,1) = 0.5*sum(temp);
end
%Shape functions for FE analysis
fem = TabShapeFnc(fem);
%Compute vectors for triplets used to assemble stiffness matrix and load
fem.ElemNDof = 2*cellfun(@length,fem.Element); % # of DOFs per element
fem.i = zeros(sum(fem.ElemNDof.^2),1); fem.j=fem.i; fem.e=fem.i; 
fem.k0 = fem.i;
fem.eDof = zeros(sum(fem.ElemNDof),1); fem.DofE = fem.eDof;
fem.iK0 = fem.i; fem.jK0 = fem.iK0; 
indexK = 0; indexF = 0; IniI = 0; IniJ = 0;
%Here, we store the triplets to assemble global matrix and load vector
for el = 1:fem.NElem  
  NDof = fem.ElemNDof(el); NKE = NDof^2;
  eDof = reshape([2*fem.Element{el}-1;2*fem.Element{el}],NDof,1); % Element DOFs
  %Triplets for stiffness matrix
  I=repmat(eDof ,1,NDof); J = I';
  fem.i(indexK+1:indexK+NKE) = I(:);
  fem.j(indexK+1:indexK+NKE) = J(:); 
  fem.eDof(indexF+1:indexF+NDof) = eDof;
  fem.e(indexK+1:indexK+NKE) = el;
  fem.DofE(indexF+1:indexF+NDof) = el; 
  %Local stiffness matrix for linear material
  Ke = LocalK(fem,fem.Element{el});
  fem.k0(indexK+1:indexK+NKE) = Ke(:);
  %Triplets to assemble element-wise stiffness matrix for solid material 
  %(used for efficient computation of stress sensitivity)
  I = repmat((1:NDof)',1,NDof); J = I';
  fem.iK0(indexK+1:indexK+NKE) = I(:)+IniI;
  fem.jK0(indexK+1:indexK+NKE) = J(:)+IniJ;
  IniI = IniI+NDof;
  IniJ = IniJ+NDof;
  indexK = indexK + NDof^2; indexF = indexF + NDof; 
end
%Compute external load vector
NLoad = size(fem.Load,1);
fem.Fext = zeros(2*fem.NNode,1);  % External load vector
if NLoad>0
  fem.Fext(2*fem.Load(1:NLoad,1)-1,1) = fem.Load(1:NLoad,2);  %x-crdnt
  fem.Fext(2*fem.Load(1:NLoad,1),1)   = fem.Load(1:NLoad,3);  %y-crdnt
end
%Obtain fixed DOFs and free DOFs
NSupp = size(fem.Supp,1);
FixedDofs = [fem.Supp(1:NSupp,2).*(2*fem.Supp(1:NSupp,1)-1);
             fem.Supp(1:NSupp,3).*(2*fem.Supp(1:NSupp,1))];
FixedDofs = FixedDofs(FixedDofs>0);
AllDofs   = 1:2*fem.NNode;
fem.FreeDofs = setdiff(AllDofs,FixedDofs);
%Compute matrices B at the centroid of each element (for VM stress computation)
idx = 0; 
B0 = zeros(sum(fem.ElemNDof)*3,1); irow = B0; icol = irow; IniI = 0; IniJ = 0;
for el = 1:fem.NElem  
  nn = fem.ElemNDof(el)/2; % Number of nodes for element e
  eNode = fem.Element{el}; % Nodes defining element e
  [~,dNdxi] = PolyShapeFnc(nn,[0,0]); % Derivative of shape function at centroid
  J0 = fem.Node(eNode,:)'*dNdxi; % Jacobian at centroid
  dNdx = dNdxi/J0;
  B = zeros(3,2*nn);
  B(1,1:2:2*nn) = dNdx(:,1)';  B(2,2:2:2*nn) = dNdx(:,2)';
  B(3,1:2:2*nn) = dNdx(:,2)';  B(3,2:2:2*nn) = dNdx(:,1)';
  B0(idx+1:idx+6*nn) = B(:);
  I = repmat(1:3,1,2*nn); J = repmat(1:2*nn,3,1);
  irow(idx+1:idx+6*nn) = I(:) + IniI;  icol(idx+1:idx+6*nn) = J(:) + IniJ;
  idx = idx + 6*nn; IniI = IniI + 3; IniJ = IniJ + 2*nn; 
end
% Assemble the Global B0 Matrix. This is a block diagonal matrix containing
% B matrices for all elements:
% B0 = [B1             ]
%      [   B2     0    ]
%      [      ...      ]
%      [   0       B_Ne]
fem.B0 = sparse(irow,icol,B0);
%Index array to assemble material tangent matrix D at element centroids
fem.rowD = zeros(9*size(fem.NElem,3),1); fem.colD = fem.rowD; 
IniI = 0; IniJ = 0; count = 1;
for el=1:fem.NElem
  I = repmat(1:3,1,3); J = repmat(1:3,3,1);
  fem.rowD(count:count+8) = I(:)+IniI;
  fem.colD(count:count+8) = J(:)+IniJ;
  IniI = IniI+3; IniJ = IniJ+3; count = count+9;
end
%Compute shape functions for FE analysis
W = cell(length(fem.ShapeFnc),1);
dNdxi = cell(length(fem.ShapeFnc),1);
for i=1:length(fem.ShapeFnc)
  if(isstruct(fem.ShapeFnc{i}))
    W{i} = fem.ShapeFnc{i}.W;
    dNdxi_temp = zeros(i,2*length(W{i}));
    for j=1:length(W{i})
      dNdxi_temp(:,(j-1)*2+1:j*2) = fem.ShapeFnc{i}.dNdxi(:,:,j);
    end
    dNdxi{i} = dNdxi_temp;
  else
    W{i} = [1];
    dNdxi{i} = [1];
  end
end
fem.W = W; fem.dNdxi = dNdxi;
%% --------------------------------------------- ELEMENT STIFFNESS MATRICES
function Ke = LocalK(fem,eNode)
E0 = fem.MatParam(1); Nu0 = E0/(2*fem.MatParam(3))-1;
D = E0/(1-Nu0^2)*[1 Nu0 0;Nu0 1 0;0 0 (1-Nu0)/2]; % Plane stress
nn = length(eNode); Ke = zeros(2*nn,2*nn); 
W = fem.ShapeFnc{nn}.W;
for q = 1:length(W)  %quadrature loop
  dNdxi = fem.ShapeFnc{nn}.dNdxi(:,:,q);
  J0 = fem.Node(eNode,:)'*dNdxi; 
  dNdx = dNdxi/J0;
  B = zeros(3,2*nn);
  B(1,1:2:2*nn) = dNdx(:,1)'; 
  B(2,2:2:2*nn) = dNdx(:,2)'; 
  B(3,1:2:2*nn) = dNdx(:,2)'; 
  B(3,2:2:2*nn) = dNdx(:,1)';
  Ke = Ke+B'*D*B*W(q)*det(J0); 
end
%% ----------------------------------------------- TABULATE SHAPE FUNCTIONS
function fem = TabShapeFnc(fem)
ElemNNode = cellfun(@length,fem.Element); % number of nodes per element
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
%% ---------------------------------------------- POLYGONAL SHAPE FUNCTIONS
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
%% -------------------------------------------------- POLYGON TRIANGULATION
function [p,Tri] = PolyTrnglt(nn,xi)
p = [cos(2*pi*((1:nn))/nn); sin(2*pi*((1:nn))/nn)]';
p = [p; xi];
Tri = zeros(nn,3); Tri(1:nn,1)=nn+1;
Tri(1:nn,2)=1:nn; Tri(1:nn,3)=2:nn+1; Tri(nn,3)=1;
%% --------------------------------------------------- POLYGONAL QUADRATURE
function [weight,point] = PolyQuad(nn)
[W,Q]= TriQuad;                  %Integration pnts & wgts for ref. triangle
[p,Tri] = PolyTrnglt(nn,[0 0]);  %Triangulate from origin
point=zeros(nn*length(W),2); weight=zeros(nn*length(W),1);
for k=1:nn
  sctr = Tri(k,:);
  for q=1:length(W)
    [N,dNds] = TriShapeFnc(Q(q,:));  %compute shape functions
    J0 = p(sctr,:)'*dNds;
    l = (k-1)*length(W) + q;
    point(l,:) = N'*p(sctr,:);
    weight(l) = det(J0)*W(q);
  end                                 
end
%% -------------------------------------------------- TRIANGULAR QUADRATURE
function [weight,point] = TriQuad
point=[1/6,1/6;2/3,1/6;1/6,2/3]; weight=[1/6,1/6,1/6];   
%% --------------------------------------------- TRIANGULAR SHAPE FUNCTIONS
function [N,dNds] = TriShapeFnc(s)
N=[1-s(1)-s(2);s(1);s(2)]; dNds=[-1,-1;1,0;0,1];
%-------------------------------------------------------------------------%