%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s1, s2, s3, s4 ]=triangulate_element(phi,coords)
    % compute the coordinates of the cut nodes
    localx=[-1 1 1 -1];
    localy=[-1 -1 1 1];
    cut_face=find(phi.*phi([2 3 4 1])<0);
    cut_val=-phi./(phi([2 3 4 1])-phi);
    x_cut=zeros(1,length(cut_face));
    y_cut=zeros(1,length(cut_face));
    for ii=1:length(cut_face)
        switch cut_face(ii)
            case 1
                y_cut(ii)=-1;
                x_cut(ii)=-1+2*cut_val(1);
            case 2
                x_cut(ii)=1;
                y_cut(ii)=-1+2*cut_val(2);
            case 3
                x_cut(ii)=1-2*cut_val(3);
                y_cut(ii)=1;
            case 4
                x_cut(ii)=-1;
                y_cut(ii)=1-2*cut_val(4);
        end
    end
    % Coordinates
    XL=[localx x_cut];
    YL=[localy y_cut];
    XG=[coords(:,1)' (x_cut+1)/2*(coords(2,1)-coords(1,1))+coords(1,1)];
    YG=[coords(:,2)' (y_cut+1)/2*(coords(3,2)-coords(2,2))+coords(2,2)];
    
    % Local levelset for nodes
    phiL=[phi zeros(size(x_cut))];
    [triangles, cellcase]=triangulate(phi); % based on marching squares

    % Generate elements
    s1 = [XG' YG'];
    s2 = triangles;
    s3 = (sum(phiL(triangles),2)>0)*1.0;
    s4 = [XL' YL'];

    % Determine ambiguity in case of cell case 5 or 10 
    if cellcase== 5 || cellcase == 10
        % disp(cellcase); 
        % Switch case if solid/void mismatch
        s3(5:6) = (sum(phi) > 0) *1.0;
    end
