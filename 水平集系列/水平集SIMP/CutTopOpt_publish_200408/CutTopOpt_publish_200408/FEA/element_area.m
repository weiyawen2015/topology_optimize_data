%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ar2=element_area(XYL,triangles,triangle_phase)
% Computes the element area contribution from each subtriangle of mother
% element
XL=XYL(:,1)';
YL=XYL(:,2)';

% Initialize arrays
phase=zeros(length(triangles)*3,1);
areaT=zeros(length(triangles),1);

% Loop over triangles to find gp, weights and areas
for ii=1:length(triangles)
    xTL = XL(triangles(ii,:));
    yTL = YL(triangles(ii,:));
    phase(ii*3-2:ii*3)=triangle_phase(ii)*[1 1 1]';
    % Triangle area
    areaT(ii) = abs(xTL(1)*yTL(2)-yTL(1)*xTL(2)+ xTL(2)*yTL(3)-yTL(2)*xTL(3)+xTL(3)*yTL(1)-yTL(3)*xTL(1))/2;
end
ar2 = areaT'*triangle_phase/4;