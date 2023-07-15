function [nodes,eles,eleNum,nodeNum] = GenerateMesh(lx,ly,numx,numy)
% To Mesh a membrane using 4 noded Elements        |
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Purpose:
%         To Mesh a square/rectangular membrane to use in FEM Analysis
% Synopsis :
%         [coordinates,nodes,nel,nnode] = GenerateMesh(lx,ly,numx,numy)
% Variable Description:
% Input :
%           lx  - length of the membrane along X-axes
%           ly  - breadth of the membrane along Y-axes
%           numx - number of Elements along X-axes
%           numy - number of Elements along Y-axes
% Output :
%           nodes - The nodal coordinates of the mesh
%           -----> nodes = [node X Y]
%           eles - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]
%           nel- total number of elements
%           nnodes- total number of nodes
%--------------------------------------------------------------------------
% lx = numx;
% ly = numy;

eleNum = numx*numy ;        % Total Number of Elements in the Mesh
nnel = 4 ;           % Number of nodes per Element
% Number of points on the Length and Breadth
npx = numx+1 ;
npy = numy+1 ;
nodeNum = npx*npy ;      % Total Number of Nodes in the Mesh
% Discretizing the Length and Breadth of the membrane
nx = linspace(0,lx,npx) ;
ny = linspace(0,ly,npy) ;
[xx,yy] = meshgrid(nx,fliplr(ny));
XX=xx;YY=yy;
nodes = [XX(:) YY(:)];

% To get the Nodal Connectivity Matrix
NodeNo = 1:nodeNum ;
eles = zeros(eleNum,nnel) ;

% small code for selecting nodal point number for specific element
NodeNo = reshape(NodeNo,npy,npx);
eles(:,4) = reshape(NodeNo(1:npy-1,1:npx-1),eleNum,1);
eles(:,3) = reshape(NodeNo(1:npy-1,2:npx),eleNum,1);
eles(:,2) = reshape(NodeNo(2:npy,2:npx),eleNum,1);
eles(:,1) = reshape(NodeNo(2:npy,1:npx-1),eleNum,1);

%
% Plotting the Finite Element Mesh
% Initialization of the required matrices
X = zeros(nnel,eleNum) ;
Y = zeros(nnel,eleNum) ;
% Extract X,Y coordinates for the (iel)-th element
for iel = 1:eleNum
    X(:,iel) = nodes(eles(iel,:),1) ;
    Y(:,iel) = nodes(eles(iel,:),2) ;
end

nodeNum = size(nodes,1) ;
eleNum = size(eles,1) ;

% nodes = [0,1;
%          0,0.5;
%          0,0;
%          0.5,1;
%          0.5,0.5;
%          0.5,0;
%          1,1;
%          1,0.5;
%          1,0;];
%      
% eles = [2 5 4 1;
%         3 6 5 2;
%         5 8 7 4;
%         6 9 8 5;];
