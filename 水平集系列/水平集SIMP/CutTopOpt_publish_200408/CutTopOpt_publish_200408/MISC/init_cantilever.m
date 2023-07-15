%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [design,mesh,fem,opt] = init_cantilever()
% this function initializes the design variables, boundary conditions and
% external load for the cantilever beam problem

%% BACKGROUND MESH
fem.lx = 2;
fem.ly = 1;
fem.nelx = 100;
fem.nely = 50;

%% Optimization parameters
opt.volfrac=0.4;
opt.max_iter = 1000;
opt.Rnodes = 5; % Smoothing in nodes
opt.deta = 0.01;

opt.initDesign = '6by4';
% opt.initDesign = 'solid';
% opt.initDesign = 'top88';

%% Material
fem.E = 1;
fem.Emin=1e-9;
fem.nu = 0.3;

%% Asymptotes and move limit
opt.move=0.02; % Move limit

opt.asyminit=0.5;
opt.asymdecrease=0.7;
opt.asymincrease=1.2;

%% Create mesh
mesh = quadmesh(fem.nelx,fem.nely,fem.lx,fem.ly);
mesh.name = 'MBB';

%% INITIAL DESIGN VARIABLES
design = ones(mesh.nnodes,1); % all nodes outside rectangle are initialized with -0.1
% Choose initial config
switch opt.initDesign
    case '6by4'
        % Make holes
        for e = 1:mesh.nel
            edof = mesh.IX(e,:);
            for n = 1:4
                x = mesh.XY(edof(n),1); % parent element x coordinates
                y = mesh.XY(edof(n),2); % parent element y coordinates
                VAL = cos((6.0*pi*x)/fem.lx)*cos(4.0*pi*y) + 0.15;
                design(edof(n)) = VAL;
            end
        end
        design = design/2 + 0.5;
        design(design>1) = 1;
        
    case 'solid'
        design=design*0+0.55;
        
    case 'top88'
        % Map to nodes
        design = zeros((mesh.nely+1)*(mesh.nelx+1),2);
        temp=top88_model(mesh.nelx,mesh.nely,max(opt.volfrac,0.3),3,1.4,1,'cant',100,0);
        temp=flipud(temp);
        elem=0;
        for ex=1:mesh.nelx
            for ey=1:mesh.nely
                elem=elem+1;
                n1 = ey + (mesh.nely+1)*(ex-1);
                n2 = ey + (mesh.nely+1)*ex;
                ndof=[n1 n2 n2+1 n1+1];
                for i=1:4
                    design(ndof(i),1) = design(ndof(i),1) + temp(ey,ex);
                    design(ndof(i),2) = design(ndof(i),2) + 1;
                end
            end
        end
        design = design(:,1)./design(:,2);
        design = design(:);
        
    otherwise
        error('opt.initDesign not recognized\n');
end

%% BC & LOAD
clamped_nodes = find(mesh.XY(:,1)==0);
bc.nodal.type = 'nodal';
bc.nodal.xnodes = clamped_nodes;
bc.nodal.ynodes = clamped_nodes;

ff.pointload.coords = [2,0.5];
ff.pointload.direction = 'y';
ff.pointload.value = -1;

% save in fem
fem.bc = bc;
fem.ff = ff;

end




