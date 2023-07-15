%----------------------------- PolyStress --------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load] = Mesh_Crack_Prob(Ne_ap)
L = 2; % Domain size
nn = floor(round(sqrt(Ne_ap/2))); he=L/(2*nn); 
NElem = round(2*nn^2); 
[X,Y] = meshgrid(he/2:he:L/2-he/2,-L/2+he/2:he:L/2-he/2);
P = [X(:) Y(:)]; % Mesh seed
[Node,Element,Supp,Load,~] = PolyMesher(@CrackDomain,NElem,0,P);
%-------------------------------------------------------------------------%