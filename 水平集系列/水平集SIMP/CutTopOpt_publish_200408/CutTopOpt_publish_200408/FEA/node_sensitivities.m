%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dc,dv] = node_sensitivities(fem,opt,mesh,phi,rho,U)
% this function calculates the nodal design sensitivities using a
% semi-analytical method

% Unwrap input
volfrac=opt.volfrac;

% Initialize 
idx=zeros(4,length(mesh.cut_ele));
dci=zeros(4,length(mesh.cut_ele));
dvi=zeros(4,length(mesh.cut_ele));

for i=1:length(mesh.cut_ele)
    % element nodes and degrees of freedom
    e = mesh.cut_ele(i);
    enodes = mesh.IX(e,:);
    edof = reshape([enodes*2-1; enodes*2],1,8);
    Ue1 = U(edof,1);
    Ue2 = Ue1;

    xR = mesh.XY(mesh.IX(e,:),1);
    yR = mesh.XY(mesh.IX(e,:),2);
    
    % stiffness matrix - no pertubation
    XYL=mesh.subtriangles{i,4};
    triangles=mesh.subtriangles{i,2};
    triangle_phase=mesh.subtriangles{i,3};

    % compute reference element matrix
    [ke_org,~]=single_element_integration(fem,xR,yR,XYL,triangles,triangle_phase);

    % loop over each node in the element and perturb it
    dce=zeros(4,1);
    dve=zeros(4,1);
    for j = 1:4
        phi_e = phi(enodes); % original level set values
        pertb = phi_e(j)/1e6; 
        phi_e(j) = phi_e(j) + pertb; % perturb level set values
        if abs(pertb) < 1e-14
            disp('small pert in nodesens')
        end
        
        % compute perturbed stiffness matrix
        [~,~,~,s4]=triangulate_element(phi_e',[xR yR]); % only use new node positions
        [ke_pertb,rho_pertb]=single_element_integration(fem,xR,yR,s4,triangles,triangle_phase);
        
        % perform FD approx
        keds = (ke_pertb-ke_org)/pertb; % differentiated ke with finite difference
        dce(j) = - Ue1'*keds*Ue2; % sensitivities of objective
        
        fvol_org = sum(rho)/(volfrac*mesh.nel)-1;
        fvol_pertb = (sum(rho)+rho_pertb-rho(e))/(volfrac*mesh.nel)-1;
        dve(j)=  (fvol_pertb-fvol_org)/pertb; % differentiated volume constraint
    end
    % Collect contribution
    idx(:,i)=enodes';
    dci(:,i)=dce;
    dvi(:,i)=dve;
end
% Form gradients
dc=full(sparse(idx(:),ones(size(idx(:),1),1),dci(:),mesh.nnodes,1));
dv=full(sparse(idx(:),ones(size(idx(:),1),1),dvi(:),mesh.nnodes,1));















