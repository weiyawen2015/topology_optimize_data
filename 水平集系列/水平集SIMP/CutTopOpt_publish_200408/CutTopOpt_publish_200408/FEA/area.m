%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g,varargout] = area(mesh,phi)
% this function computes the areas of solid domain (and possibly sens) 
rho = zeros(mesh.nel,1);
rho(mesh.void_ele) = 0;
rho(mesh.solid_ele) = 1;

% loop over each cut element and perform triangulation
for i=1:length(mesh.cut_ele)
    e = mesh.cut_ele(i);
    %
    ar2=element_area(mesh.subtriangles{i,4},mesh.subtriangles{i,2},mesh.subtriangles{i,3});
    
    rho(e) = rho(e) + ar2;
end
g=sum(rho)/mesh.nel;

if nargout > 1
% this function calculates the nodal design sensitivities using a
% semi-analytical method

% Unwrap input

% Initialize 
dvi=zeros(4,length(mesh.cut_ele));
idx=zeros(4,length(mesh.cut_ele));

for i=1:length(mesh.cut_ele)
    % element nodes and degrees of freedom
    e = mesh.cut_ele(i);
    enodes = mesh.IX(e,:);
    xR = mesh.XY(mesh.IX(e,:),1);
    yR = mesh.XY(mesh.IX(e,:),2);
    % loop over each node in the element and perturb it
    dve=zeros(4,1);
    for j = 1:4
        phi_e = phi(enodes); % original level set values
        pertb = phi_e(j)/1e6; % small pertubation
        phi_e(j) = phi_e(j) + pertb; % perturb level set values
    
        % perturbed stiffness matrix
        [~,~,~,s4]=triangulate_element(phi_e',[xR yR]); % only use new node positions
        rho_pertb=element_area(s4,mesh.subtriangles{i,2},mesh.subtriangles{i,3});

        fvol_org = sum(rho); % total volume
        fvol_pertb = (sum(rho)+rho_pertb-rho(e)); 
        dve(j)=  (fvol_pertb-fvol_org)/pertb; % differentiated volume constraint
    end
    % Collect contribution
    idx(:,i)=enodes';
    dvi(:,i)=dve;
end
% Form gradients based on sparse assembly
varargout{1}=full(sparse(idx(:),ones(size(idx(:),1),1),dvi(:),mesh.nnodes,1))/mesh.nel;

end














