%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [k, f_NL] = Global_K_compile(ElemNDof,Element,MatModel,...
                                      MatParam,Thickness,E,W,dNdxi,Node,d)
k = zeros(sum(ElemNDof.^2),1);
f_NL = zeros(sum(ElemNDof),1); 
indexK = 0; indexF = 0;
NElem = length(E);
for el = 1:NElem
  NDof = ElemNDof(el);
  eDof = reshape([2*Element{el}-1;2*Element{el}],NDof,1); % Element DOFs
  ue = d(eDof); % Element displacement vector
  [Ke, Fint_e] = ...
    LocalK(Node,W,dNdxi,Thickness,NDof,ue,MatModel,MatParam,Element{el});
  k(indexK+1:indexK+NDof^2) = Ke(:); % Stiffness matrix of solid elements
  f_NL(indexF+1:indexF+NDof) = Fint_e(:); % For load vector
  indexK = indexK + NDof^2; indexF = indexF + NDof;
end
%% --------------------- ELEMENT STIFFNESS MATRIX AND INTERNAL FORCE VECTOR
function [Ke,Fint_e] = LocalK(Node,W,dNdxi,Thickness,NDof,ue,MatModel,...
    MatParam,eNode)
Ke = zeros(NDof,NDof); Fint_e = zeros(NDof,1); 
nn = length(eNode);
W = W{nn}; dNdxi = dNdxi{nn};
B = zeros(3,NDof); 
Node_el = Node(eNode,:)'; % Nodal coordinates of undeformed element
for q = 1:length(W) % Loop over Gauss points
  dNdxi_temp = dNdxi(:,2*(q-1)+1:2*q);
  J0 = Node_el*dNdxi_temp; WdJ = W(q)*det(J0);
  dNdX = dNdxi_temp/J0;
  B(1,1:2:NDof) = dNdX(:,1); B(2,2:2:NDof) = dNdX(:,2); 
  B(3,1:2:NDof) = dNdX(:,2); B(3,2:2:NDof) = dNdX(:,1);
  E = B*ue; % Infinitesimal strain vector
  [S, D_GP] = material_model(MatModel,MatParam,E);
  Fint_e = Fint_e + B'*S.*WdJ; % Internal force vector
  Ke = Ke + B'*D_GP*B*WdJ; % Tangent stiffness matrix
end
Ke = Ke.*Thickness; Fint_e = Fint_e.*Thickness; 
%-------------------------------------------------------------------------%