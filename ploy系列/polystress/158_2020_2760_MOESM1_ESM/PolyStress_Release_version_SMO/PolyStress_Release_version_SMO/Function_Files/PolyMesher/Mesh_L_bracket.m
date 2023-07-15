%----------------------------- PolyStress --------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load] = Mesh_L_bracket(Ne_ap)
L = 1; t = 2*L/5;
nn = 2*floor(round(sqrt(Ne_ap/4))/2); he = t/nn; 
NElem = round(4*nn^2); 
[X1,Y1] = meshgrid(he/2:he:L-he/2,he/2:he:t-he/2); 
[Y2,X2] = meshgrid(t+he/2:he:L-he/2,he/2:he:t-he/2);
P = [X1(:) Y1(:);X2(:) Y2(:)]; % Mesh seed
[Node,Element,Supp,Load,~] = PolyMesher(@L_bracketDomain,NElem,0,P);
%-------------------------------------------------------------------------%