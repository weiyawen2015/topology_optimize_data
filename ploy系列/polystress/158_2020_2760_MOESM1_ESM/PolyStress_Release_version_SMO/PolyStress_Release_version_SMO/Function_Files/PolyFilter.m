%-------------------------------- PolyTop --------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyTop: A Matlab %
% implementation of a general topology optimization framework using       %
% unstructured polygonal finite element meshes", Struct Multidisc Optim,  %
% DOI 10.1007/s00158-011-0696-x                                           %
%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [P] = PolyFilter(fem,R,q,Sym)
if R<=0, P = speye(fem.NElem); return; end %P is set to identity when R<=0
ElemCtrd = zeros(fem.NElem,2);
for el = 1:fem.NElem        %Compute the centroids of all the elements
  vx = fem.Node(fem.Element{el},1); vy = fem.Node(fem.Element{el},2);
  temp = vx.*vy([2:end 1])-vy.*vx([2:end 1]);
  A = 0.5*sum(temp);
  ElemCtrd(el,1) = 1/(6*A)*sum((vx+vx([2:end 1])).*temp);
  ElemCtrd(el,2) = 1/(6*A)*sum((vy+vy([2:end 1])).*temp);
end
if nargin==4
  switch Sym
    case 'X',  x = [ElemCtrd(:,1),abs(ElemCtrd(:,2))];
    case 'Y',  x = [abs(ElemCtrd(:,1)),ElemCtrd(:,2)];
    case 'XY', x = [abs(ElemCtrd(:,1)),abs(ElemCtrd(:,2))];            
    otherwise, error('Use Sym = "X", Sym = "Y", or Sym = "XY"');
  end
else
  x = ElemCtrd;
end
[d] = DistPntSets(x,x,R); %Obtain distance values & indexes
P = sparse(d(:,1),d(:,2),(1-d(:,3)/R).^q);   %Assemble filter matrix
P = spdiags(1./sum(P,2),0,fem.NElem,fem.NElem)*P;
%---------------------------------- COMPUTE DISTANCE BETWEEN TWO POINT SETS
function [d] = DistPntSets(PS1,PS2,R)
d = cell(size(PS1,1),1);
for el = 1:size(PS1,1)       %Compute the distance information
  dist = sqrt((PS1(el,1)-PS2(:,1)).^2 + (PS1(el,2)-PS2(:,2)).^2);
  [I,J] = find(dist<=R);   %Find the indices for distances less that R
  d{el} = [I,J+(el-1),dist(I)];
end
d = cell2mat(d);             %Matrix of indices and distance value
%-------------------------------------------------------------------------%