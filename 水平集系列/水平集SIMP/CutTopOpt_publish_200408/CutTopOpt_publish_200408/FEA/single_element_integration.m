%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make subtriangulation of a single element and create stiffness matrix
% Returns stiffness matrix and area-ratio
function [ke,areaRatio]=single_element_integration(fem,xR,yR,XYL,triangles,triangle_phase)
% Material
E=fem.E;
Emin=fem.Emin;
nu=fem.nu;

XL=XYL(:,1)';
YL=XYL(:,2)';

% Initialize arrays
xiR=zeros(size(triangles,1)*3,1);
etaR=zeros(size(triangles,1)*3,1);
wR=zeros(size(triangles,1)*3,1);
phase=zeros(size(triangles,1)*3,1);
areaT=zeros(size(triangles,1),1);

% Loop over triangles to find gp, weights and areas
for ii=1:size(triangles,1)
    xTL = XL(triangles(ii,:));
    yTL = YL(triangles(ii,:));
    % Find the Gauss points
    [xiR(ii*3-2:ii*3),etaR(ii*3-2:ii*3),wR(ii*3-2:ii*3)] = subtriangle_gaussPoints2(xTL,yTL); % dep on levelset
    % Phase
    phase(ii*3-2:ii*3)=triangle_phase(ii)*[1 1 1]';
    % Triangle area
    areaT(ii) = abs(xTL(1)*yTL(2)-yTL(1)*xTL(2)+ xTL(2)*yTL(3)-yTL(2)*xTL(3)+xTL(3)*yTL(1)-yTL(3)*xTL(1))/2;
    % Triangle weight
    wR(ii*3-2:ii*3)=wR(ii*3-2:ii*3)*areaT(ii); %
end
areaRatio = areaT'*triangle_phase/4;

% Constant quantities for integration
Cmat=1/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]; % constitutive matrix for plane stress
L = zeros(3,4);
L(1,1) = 1; L(2,4) = 1; L(3,2:3) = 1; % strain mapping vector

% Perform integration over rectangular element using Gauss quadrature
ke = zeros(8,8);
for i = 1:length(xiR)
    xi = xiR(i);
    eta = etaR(i);
    [~,dNxi,dNeta,~,dNmat] = shapeFunc2D('Q1',xi,eta); % dep on GP
    J = [dNxi; dNeta]*[xR,yR]; % jacobian
    detJ = det(J);
    invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
    G = [invJ zeros(2); zeros(2) invJ]; % mapping matrix between reference element and global system
    B = L*G*dNmat; % strain-displacement matrix
    weight = wR(i)*detJ; % gauss weight
    ke = ke + (Emin+(E-Emin)*phase(i))*(B'*Cmat*B)*weight;
end

end

function [xiR,etaR,wR] = subtriangle_gaussPoints2(xT,yT)
% let the xT and yT be given in the rectangular -1;1 x -1;1
xiT  = [1/6 2/3 1/6];
etaT = [1/6 1/6 2/3];
wT  = [1/3 1/3 1/3];

xiR=zeros(size(xiT));
etaR=zeros(size(xiT));
% compute gauss points in reference rectangle coordinates
for i = 1:length(xiT)
    N = [1-xiT(i)-etaT(i), xiT(i), etaT(i)];
    xiR(i) = xT(1)*N(1)+xT(2)*N(2)+xT(3)*N(3);
    etaR(i) = yT(1)*N(1)+yT(2)*N(2)+yT(3)*N(3);
end
wR = wT; % weights remain the same

end

