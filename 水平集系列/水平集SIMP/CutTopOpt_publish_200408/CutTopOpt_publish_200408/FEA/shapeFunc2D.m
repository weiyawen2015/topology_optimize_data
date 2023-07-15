%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,dNxi,dNeta,Nmat,dNmat] = shapeFunc2D(eleType,xi,eta)
 % this function returns the shape function and differentiated shape
% function matrix for a triangular or quadrilateral element

if strcmp(eleType,'Q1')
    % Shape functions
    N = 1/4*[(1-eta)*(1-xi),(1-eta)*(1+xi),(1+eta)*(1+xi),(1-xi)*(1+eta)];
    % Differentiated shape functions
    dNxi = 1/4*[-(1-eta)  (1-eta) (1+eta) -(1+eta)];
    dNeta = 1/4*[-(1-xi) -(1+xi) (1+xi)  (1-xi)];
    
    % Shape functions matrices
    Nmat = [N(1), 0, N(2), 0, N(3), 0, N(4) 0;
        0, N(1), 0, N(2), 0, N(3), 0, N(4)];
    
    % differentiate shape functions - gradients   
    dNmat = 1/4*[-(1-eta) 0 (1-eta) 0 (1+eta) 0 -(1+eta) 0
                  -(1-xi) 0 -(1+xi) 0 (1+xi) 0 (1-xi) 0
                 0 -(1-eta) 0 (1-eta) 0 (1+eta) 0 -(1+eta)
                 0 -(1-xi) 0 -(1+xi) 0 (1+xi) 0 (1-xi)];
else
    disp('Not implemented');
end

end

