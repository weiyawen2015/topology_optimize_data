%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = quadmesh(nelx,nely,lx,ly)
% this function creates a quadrilateral background mesh and returns the
% associated node positions, facets and element/node connectivity tables

dx=lx/nelx;
dy=ly/nely;
[x,y] = meshgrid(0:dx:lx,0:dy:ly);
XY = [x(:) y(:)];
nel = nelx*nely;
nnodes = size(XY,1);
IX = eleConnectivity(nelx,nely);
nIX = nodeConnectivity(IX,nnodes);

mesh.XY = XY; % node positions
mesh.IX = IX; % element connectivity
mesh.nIX = nIX; % node connectivity
mesh.nelx = nelx;
mesh.nely = nely;
mesh.dx = dx;
mesh.dy = dy;
mesh.nel = nel;
mesh.nnodes = nnodes;

end

function IX = eleConnectivity(nelx,nely)

IX = zeros(nelx*nely,4);
row=0;
for i=1:nelx
    for j=1:nely
        row = row + 1;
        IX(row,:)=[
            (i-1)*(nely+1)+(j) ...
            (i+0)*(nely+1)+(j) ...
            (i+0)*(nely+1)+(j+1)...
            (i-1)*(nely+1)+(j+1)...
            ];
    end
end

end

function nIX = nodeConnectivity(IX,nnodes)

nIX = cell(nnodes,2);
for node=1:nnodes
    [eleIDs,~]=find(node==IX); % find the elements connected to the node
    nIX{node,1}=length(eleIDs); % number of elements connected to the node
    nIX{node,2}=eleIDs; % element id's of the connected elements
end

end

