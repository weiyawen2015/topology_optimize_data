%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function updates the mesh structure to contain information regarding
% cut elements, intersections, coordinates and connectivity
function mesh = categorize_elements(mesh,phi)
mesh.subtriangles={}; % clear the list to avoid leftovers

% Provides the categorization of elements based on the levelset phi
mesh.cut_faces=phi(mesh.IX).*phi(mesh.IX(:,[2 3 4 1]))<0; % Indication of cut faces
mesh.cut_ele=find(sum(mesh.cut_faces,2)>1); % contains all elements that are cut
mesh.solid_ele=find(sum(mesh.cut_faces,2)==0 & phi(mesh.IX(:,1))>0);
mesh.void_ele=find(sum(mesh.cut_faces,2)==0 & phi(mesh.IX(:,1))<0);
mesh.dcut_ele=find(sum(mesh.cut_faces,2)>2); 

% Compute all intersections of faces
mesh.cut_val=-phi(mesh.IX)./(phi(mesh.IX(:,[2 3 4 1]))-phi(mesh.IX));

% For all cut elements; compute triangulation
for kk=1:length(mesh.cut_ele)
    elem=mesh.cut_ele(kk);
    coords=mesh.XY(mesh.IX(elem,:),:);
    
    % Perform triangulation
    [s1,s2,s3,s4]=triangulate_element(phi(mesh.IX(elem,:))',coords);
    
    % Add to mesh struct
    mesh.subtriangles{kk,1}=s1; % Global coords
    mesh.subtriangles{kk,2}=s2; % Connectivity
    mesh.subtriangles{kk,3}=s3; % material Phase
    mesh.subtriangles{kk,4}=s4; % Rectangular ref. coordinates
end
   