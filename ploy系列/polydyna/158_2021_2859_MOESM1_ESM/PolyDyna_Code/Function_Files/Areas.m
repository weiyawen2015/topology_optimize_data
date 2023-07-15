%-------------------------------- Areas ----------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    %
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
function ElemArea = Areas(fem)
ElemArea = zeros(fem.NElem,1); 
for el = 1:fem.NElem 
  vx = fem.Node(fem.Element{el},1); vy = fem.Node(fem.Element{el},2);
  temp = vx.*vy([2:end 1])-vy.*vx([2:end 1]);
  ElemArea(el,1) = 0.5*sum(temp);
end
%-------------------------------------------------------------------------%