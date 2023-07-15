%------------------------------ Centroids --------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    %
% Structural and Multidisciplinary Optimization,                          %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
function ElemCtrd = Centroids(fem)
ElemCtrd = zeros(fem.NElem,2); 
for el = 1:fem.NElem        %Compute the centroids of all the elements
  vx=fem.Node(fem.Element{el},1); vy=fem.Node(fem.Element{el},2);
  temp = vx.*vy([2:end 1])-vy.*vx([2:end 1]);
  A = 0.5*sum(temp);
  ElemCtrd(el,1) = 1/(6*A)*sum((vx+vx([2:end 1])).*temp);
  ElemCtrd(el,2) = 1/(6*A)*sum((vy+vy([2:end 1])).*temp);
end
end
%-------------------------------------------------------------------------%