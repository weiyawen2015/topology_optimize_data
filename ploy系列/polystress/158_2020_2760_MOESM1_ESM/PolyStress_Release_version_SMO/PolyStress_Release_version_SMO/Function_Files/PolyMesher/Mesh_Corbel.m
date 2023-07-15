%----------------------------- PolyStress --------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load] = Mesh_Corbel(Ne_ap)
L = 2; % Domain size
nn = floor(round(0.5*sqrt(Ne_ap))); he=L/nn;
NElem = round(4*nn^2); 
[X1,Y1] = meshgrid(he/2:he:L-he/2,-1.5*L+he/2:he:1.5*L-he/2);
[X2,Y2] = meshgrid(L+he/2:he:2*L-he/2,-0.5*L+he/2:he:0.5*L-he/2);
P = [X1(:) Y1(:);X2(:) Y2(:)]; % Mesh seed
[Node,Element,Supp,Load,~] = PolyMesher(@CorbelDomain,NElem,0,P);
%-------------------------------------------------------------------------%