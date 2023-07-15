%------------------------------ PolyDyna ---------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    % 
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load] = Mesh_Support(Ne_ap)
L = 3;
Ny = floor(round(sqrt(Ne_ap)));
Nx = round(Ny); if mod(Nx,2)~=0, Nx = Nx+1; end
hx = L/Nx; % Estimated element size in x direction
hy = L/Ny; % Estimated element size in y direction
NElem = round(Nx*Ny); % Number of elements
[X,Y] = meshgrid(-L/2+hx/2:hx:L/2-hx/2,hy/2:hy:L-hy/2); 
P = [X(:) Y(:)]; % Mesh seed
Domain = @SupportDomain;
[Node,Element,Supp,Load,~] = PolyMesher(Domain,NElem,0,P);
%-------------------------------------------------------------------------%