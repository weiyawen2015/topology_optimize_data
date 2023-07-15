function NURBS = Geom_Mod(L, W, Order, Num, BoundCon)
switch BoundCon
    case {1, 2, 3}
        knots{1} = [0 0 1 1]; knots{2} = [0 0 1 1];
        ControlPts(:,:,1) = [0 L; 0 0; 0 0; 1 1];
        ControlPts(:,:,2) = [0 L; W W; 0 0; 1 1];
    case 4
        knots{1} = [0 0 0.5 1 1]; knots{2} = [0 0 1 1];
        ControlPts(:,:,1) = [0 0 L; L 0 0; 0 0 0; 1 1 1];
        ControlPts(:,:,2) = [W W L; L W W; 0 0 0; 1 1 1];
    case 5
        W = W/2;
        knots{1} = [0 0 0 1 1 1]; knots{2} = [0 0 1 1];
        ControlPts(:,:,1) = [0 W W; W W 0; 0 0 0; 1 sqrt(2)/2 1];
        ControlPts(:,:,2) = [0 L L; L L 0; 0 0 0; 1 sqrt(2)/2 1];
end
coefs = zeros(size(ControlPts));
coefs(1,:,:) = ControlPts(1,:,:).*ControlPts(4,:,:);
coefs(2,:,:) = ControlPts(2,:,:).*ControlPts(4,:,:);
coefs(3,:,:) = ControlPts(3,:,:).*ControlPts(4,:,:);
coefs(4,:,:) = ControlPts(4,:,:);
NURBS = nrbmak(coefs, knots);
NURBS = nrbdegelev(NURBS,Order);
nrbplot(NURBS,[100 100],'light','on')
iknot_u = linspace(0,1,Num(1)); iknot_v = linspace(0,1,Num(2));
NURBS = nrbkntins(NURBS,{setdiff(iknot_u,NURBS.knots{1}),setdiff(iknot_v,NURBS.knots{2})});
end
%======================================================================================================================%
% Subfunction Geom_Mod:                                                                                                %
%                                                                                                                      %
% A compact and efficient MATLAB code for the numerical implementation isogeometric topology optimization in 2D        %
%                                                                                                                      %
% Developed by: Jie Gao                                                                                                %
% Email: JieGao@hust.edu.cn                                                                                            %
%                                                                                                                      %
% Main references:                                                                                                     %
%                                                                                                                      %
% (1) Jie Gao, Lin Wang, Zhen Luo, Liang Gao. IgaTop: an implementation of topology optimization for structures        %
% using IGA in Matlab. Structural and multidisciplinary optimization.                                                  %
%                                                                                                                      %
% (2) Jie Gao, Liang Gao, Zhen Luo, Peigen Li. Isogeometric topology optimization for continuum structures using       %
% density distribution function. Int J Numer Methods Eng, 2019, 119:991¨C1017                                          %
%                                                                                                                      %
% *********************************************   Disclaimer   ******************************************************* %
% The authors reserve all rights for the programs. The programs may be distributed and used for academic and           %
% educational purposes. The authors do not guarantee that the code is free from errors,and they shall not be liable    %
% in any event caused by the use of the program.                                                                       %
%======================================================================================================================%