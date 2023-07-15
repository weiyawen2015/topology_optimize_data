function [DBoudary, F] = Boun_Cond(CtrPts, BoundCon, NURBS, Dofs_Num)
%% boundary conditions
switch BoundCon
    case 1 % Cantilever beam
        DBoudary.CtrPtsOrd = CtrPts.Seque(:,1);            % Dirichlet boundary conditions
        load.u = 1; load.v = 0.5;
        [N, id] = nrbbasisfun([load.u; load.v], NURBS);
        NBoudary.CtrPtsOrd = id'; NBoudary.N = N;          % Neumann boundary conditions
    case 2 % MBB beam
        DBoudary.CtrPtsOrd1 = CtrPts.Seque(1,1); DBoudary.CtrPtsOrd2 = CtrPts.Seque(1,end);
        load.u = 0.5; load.v = 1;
        [N, id] = nrbbasisfun([load.u; load.v], NURBS);
        NBoudary.CtrPtsOrd = id'; NBoudary.N = N;
    case 3 % Michell-type structure
        DBoudary.CtrPtsOrd1 = CtrPts.Seque(1,1); DBoudary.CtrPtsOrd2 = CtrPts.Seque(1,end);
        load.u = 0.5; load.v = 0;
        [N, id] = nrbbasisfun([load.u; load.v], NURBS);
        NBoudary.CtrPtsOrd = id'; NBoudary.N = N;
    case 4 % L beam
        DBoudary.CtrPtsOrd = CtrPts.Seque(:,1);
        load.u = 1; load.v = 1;
        [N, id] = nrbbasisfun([load.u; load.v], NURBS);
        NBoudary.CtrPtsOrd = id'; NBoudary.N = N;
    case 5 % A quarter annulus
        DBoudary.CtrPtsOrd = CtrPts.Seque(:,end);
        load.u = 0; load.v = 1;
        [N, id] = nrbbasisfun([load.u; load.v], NURBS);
        NBoudary.CtrPtsOrd = id'; NBoudary.N = N;
end
%% the force imposed on the structure
F = zeros(Dofs_Num,1);
switch BoundCon
    case {1,2,3,4}
        F(NBoudary.CtrPtsOrd+CtrPts.Num) = -1*NBoudary.N;
    case 5
        F(NBoudary.CtrPtsOrd) = -1*NBoudary.N;
end
end
%======================================================================================================================%
% Subfunction Boun_Cond:                                                                                               %
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