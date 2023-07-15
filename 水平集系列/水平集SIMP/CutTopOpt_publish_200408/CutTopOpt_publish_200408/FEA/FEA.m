%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,U,rho,mesh,sedens] = FEA(fem,mesh)
% this function drives the finite element analysis by assembling the
% global stiffness matrix, and applying ghost penalties, boundary conditions
% and external loads
E=fem.E;
Emin=fem.Emin;
nu=fem.nu;

%% ASSEMBLY - UNCUT ELEMENTS
a = mesh.dx;
b = mesh.dy;
ke =  1/(1-nu^2) * [(((-nu + 1) * a ^ 2 + 2 * b ^ 2) / a / b) / 0.6e1 0.1e1 / 0.8e1 + nu / 0.8e1 (((-nu + 1) * a ^ 2 - 4 * b ^ 2) / a / b) / 0.12e2 -0.1e1 / 0.8e1 + 0.3e1 / 0.8e1 * nu (((-1 + nu) * a ^ 2 - 2 * b ^ 2) / a / b) / 0.12e2 -0.1e1 / 0.8e1 - nu / 0.8e1 (((-1 + nu) * a ^ 2 + b ^ 2) / a / b) / 0.6e1 0.1e1 / 0.8e1 - 0.3e1 / 0.8e1 * nu; 0.1e1 / 0.8e1 + nu / 0.8e1 (((-nu + 1) * b ^ 2 + 2 * a ^ 2) / a / b) / 0.6e1 0.1e1 / 0.8e1 - 0.3e1 / 0.8e1 * nu (((-1 + nu) * b ^ 2 + a ^ 2) / a / b) / 0.6e1 -0.1e1 / 0.8e1 - nu / 0.8e1 (((-1 + nu) * b ^ 2 - 2 * a ^ 2) / a / b) / 0.12e2 -0.1e1 / 0.8e1 + 0.3e1 / 0.8e1 * nu (((-nu + 1) * b ^ 2 - 4 * a ^ 2) / a / b) / 0.12e2; (((-nu + 1) * a ^ 2 - 4 * b ^ 2) / a / b) / 0.12e2 0.1e1 / 0.8e1 - 0.3e1 / 0.8e1 * nu (((-nu + 1) * a ^ 2 + 2 * b ^ 2) / a / b) / 0.6e1 -0.1e1 / 0.8e1 - nu / 0.8e1 (((-1 + nu) * a ^ 2 + b ^ 2) / a / b) / 0.6e1 -0.1e1 / 0.8e1 + 0.3e1 / 0.8e1 * nu (((-1 + nu) * a ^ 2 - 2 * b ^ 2) / a / b) / 0.12e2 0.1e1 / 0.8e1 + nu / 0.8e1; -0.1e1 / 0.8e1 + 0.3e1 / 0.8e1 * nu (((-1 + nu) * b ^ 2 + a ^ 2) / a / b) / 0.6e1 -0.1e1 / 0.8e1 - nu / 0.8e1 (((-nu + 1) * b ^ 2 + 2 * a ^ 2) / a / b) / 0.6e1 0.1e1 / 0.8e1 - 0.3e1 / 0.8e1 * nu (((-nu + 1) * b ^ 2 - 4 * a ^ 2) / a / b) / 0.12e2 0.1e1 / 0.8e1 + nu / 0.8e1 (((-1 + nu) * b ^ 2 - 2 * a ^ 2) / a / b) / 0.12e2; (((-1 + nu) * a ^ 2 - 2 * b ^ 2) / a / b) / 0.12e2 -0.1e1 / 0.8e1 - nu / 0.8e1 ...
    (((-1 + nu) * a ^ 2 + b ^ 2) / a / b) / 0.6e1 0.1e1 / 0.8e1 - 0.3e1 / 0.8e1 * nu (((-nu + 1) * a ^ 2 + 2 * b ^ 2) / a / b) / 0.6e1 0.1e1 / 0.8e1 + nu / 0.8e1 (((-nu + 1) * a ^ 2 - 4 * b ^ 2) / a / b) / 0.12e2 -0.1e1 / 0.8e1 + 0.3e1 / 0.8e1 * nu; -0.1e1 / 0.8e1 - nu / 0.8e1 (((-1 + nu) * b ^ 2 - 2 * a ^ 2) / a / b) / 0.12e2 -0.1e1 / 0.8e1 + 0.3e1 / 0.8e1 * nu (((-nu + 1) * b ^ 2 - 4 * a ^ 2) / a / b) / 0.12e2 0.1e1 / 0.8e1 + nu / 0.8e1 (((-nu + 1) * b ^ 2 + 2 * a ^ 2) / a / b) / 0.6e1 0.1e1 / 0.8e1 - 0.3e1 / 0.8e1 * nu (((-1 + nu) * b ^ 2 + a ^ 2) / a / b) / 0.6e1; (((-1 + nu) * a ^ 2 + b ^ 2) / a / b) / 0.6e1 -0.1e1 / 0.8e1 + 0.3e1 / 0.8e1 * nu (((-1 + nu) * a ^ 2 - 2 * b ^ 2) / a / b) / 0.12e2 0.1e1 / 0.8e1 + nu / 0.8e1 (((-nu + 1) * a ^ 2 - 4 * b ^ 2) / a / b) / 0.12e2 0.1e1 / 0.8e1 - 0.3e1 / 0.8e1 * nu (((-nu + 1) * a ^ 2 + 2 * b ^ 2) / a / b) / 0.6e1 -0.1e1 / 0.8e1 - nu / 0.8e1; 0.1e1 / 0.8e1 - 0.3e1 / 0.8e1 * nu (((-nu + 1) * b ^ 2 - 4 * a ^ 2) / a / b) / 0.12e2 0.1e1 / 0.8e1 + nu / 0.8e1 (((-1 + nu) * b ^ 2 - 2 * a ^ 2) / a / b) / 0.12e2 -0.1e1 / 0.8e1 + 0.3e1 / 0.8e1 * nu (((-1 + nu) * b ^ 2 + a ^ 2) / a / b) / 0.6e1 -0.1e1 / 0.8e1 - nu / 0.8e1 (((-nu + 1) * b ^ 2 + 2 * a ^ 2) / a / b) / 0.6e1;];

% sparse assembly preallocation
ndof_ele = 8;
iK = zeros(mesh.nel*ndof_ele*ndof_ele,1);
jK = zeros(mesh.nel*ndof_ele*ndof_ele,1);
KE = zeros(mesh.nel*ndof_ele*ndof_ele,1);
nK = 0;

rho = zeros(mesh.nel,1);
rho(mesh.void_ele) = 0;
rho(mesh.solid_ele) = 1;

% loop over all other elements than the cut ones
for e=setdiff(1:mesh.nel,mesh.cut_ele)
    enodes = mesh.IX(e,:);
    edof = zeros(1,8);
    edof(1:2:end) = enodes*2-1;
    edof(2:2:end) = enodes*2;
    for krow = 1:size(ke,1)
        for kcol = 1:size(ke,2)
            nK = nK+1;
            iK(nK) = edof(krow);
            jK(nK) = edof(kcol);
            KE(nK) = (Emin+(E-Emin)*rho(e))*ke(krow,kcol);
        end
    end
end

%% ASSEMBLY - CUT ELEMENTS
% loop over each cut element and perform triangulation
kee=cell(1,length(mesh.cut_ele));
for i=1:length(mesh.cut_ele)
    e = mesh.cut_ele(i);
    enodes = mesh.IX(e,:);
    edof = zeros(1,8);
    edof(1:2:end) = enodes*2-1;
    edof(2:2:end) = enodes*2;
    xR = mesh.XY(enodes,1); % parent element x coordinates
    yR = mesh.XY(enodes,2); % parent element y coordinates
    XYL=mesh.subtriangles{i,4};
    triangles=mesh.subtriangles{i,2};
    triangle_phase=mesh.subtriangles{i,3};
    
    % perform sub element integration
    [ke2,ar2]=single_element_integration(fem,xR,yR,XYL,triangles,triangle_phase);
    
    % Add to triplets
    rho(e) = rho(e) + ar2;
    for krow = 1:size(ke2,1)
        for kcol = 1:size(ke2,2)
            nK = nK+1;
            iK(nK) = edof(krow);
            jK(nK) = edof(kcol);
            KE(nK) = ke2(krow,kcol);
        end
    end
    kee{i}=ke2; % Save stiffness matrices
end
% Assemble matrix
K = sparse(iK,jK,KE,mesh.nnodes*2,mesh.nnodes*2);
K = (K+K')/2; % Symmetrize

%% BC + LOAD
% define load vector
F = zeros(2*mesh.nnodes,1);

bc_names = fieldnames(fem.bc);
n_bc = numel(bc_names);
% loop over all applied BC's
for i = 1:n_bc
    % identify BC type from bc struct field name
    bcN = fem.bc.(bc_names{i});
    if strcmp(bcN.type,'nodal') % apply BC using zero-one method
        % xnodes
        for j=1:length(bcN.xnodes)
            node=bcN.xnodes(j);
            xdof=node*2-1;
            K(:,xdof)=0;
            K(xdof,:)=0;
            K(xdof,xdof)=1;
            F(xdof,1)=0;
        end
        % ynodes
        for j=1:length(bcN.ynodes)
            node=bcN.ynodes(j);
            ydof=node*2;
            K(:,ydof)=0;
            K(ydof,:)=0;
            K(ydof,ydof)=1;
            F(ydof,1)=0;
        end
    else
        disp('Only nodal bc supported');
    end
end

% apply point load
ff_names = fieldnames(fem.ff);
n_ff = numel(ff_names);
% loop over point loads
for i = 1:n_ff
    pointLoad = fem.ff.(ff_names{i});
    % loop over all elements and check in which element the point load
    % should be applied as the point load coords do not have to coincide
    % with existing nodes in the background mesh
    for e=1:mesh.nel
        enodes = mesh.IX(e,:);
        xe = mesh.XY(enodes,1);
        ye = mesh.XY(enodes,2);
        in_poly = inpolygon(pointLoad.coords(1),pointLoad.coords(2),xe,ye);
        if logical(in_poly) == true
            xR = xe;
            yR = ye;
            if strcmp(pointLoad.direction,'x')
                loadDofs = enodes*2-1;
            elseif strcmp(pointLoad.direction,'y')
                loadDofs = enodes*2;
            end
        end
    end
    % transform global gauss points to reference rectangle
    xiR=(pointLoad.coords(1)-min(xR))*2/(max(xR)-min(xR))-1;
    etaR=(pointLoad.coords(2)-min(yR))*2/(max(yR)-min(yR))-1;
    
    % distribute load on the four nodes in the element
    [N,~,~,~,~] = shapeFunc2D('Q1',xiR,etaR);
    F(loadDofs,1) = F(loadDofs,1) + N'*pointLoad.value;
end

mesh.F = F;

%% SOLVE SYSTEM
U = K\F;

%% Compute the strain energy density of each element
sedens=size(rho);
for e=setdiff(1:mesh.nel,mesh.cut_ele)
    enodes = mesh.IX(e,:);
    edof = zeros(1,8);
    edof(1:2:end) = enodes*2-1;
    edof(2:2:end) = enodes*2;
    sedens(e)=U(edof,1)'*(Emin+(E-Emin)*rho(e))*ke*U(edof,1);
end
for i=1:length(mesh.cut_ele)
    e = mesh.cut_ele(i);
    enodes = mesh.IX(e,:);
    edof = zeros(1,8);
    edof(1:2:end) = enodes*2-1;
    edof(2:2:end) = enodes*2;
    sedens(e)=U(edof,1)'*kee{i}*U(edof,1);
end

end























