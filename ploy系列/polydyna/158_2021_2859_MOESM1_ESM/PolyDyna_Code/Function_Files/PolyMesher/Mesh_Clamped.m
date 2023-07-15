%------------------------------ PolyDyna ---------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    % 
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load] = Mesh_Clamped(Ne_ap)
L = 12; t = 2;
Ny = floor(round(sqrt(Ne_ap*t/L)));
he = t/Ny; % Estimated element size (nn must be an even number)
Nx = Ny*L/t; if mod(Nx,2)~=0, Nx = Nx+1; end
NElem = round(Nx*Ny); % Number of elements
[X,Y] = meshgrid(-L/2+he/2:he:L/2-he/2,he/2:he:t-he/2); 
P = [X(:) Y(:)]; % Mesh seed
Domain = @ClampedDomain;
[Node,Element,Supp,Load,~] = PolyMesher(Domain,NElem,0,P);
%-------------------------------------------------------------------------%