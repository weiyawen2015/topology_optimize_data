function U = Solving(CtrPts, DBoudary, Dofs, K, F, BoundCon)
switch BoundCon
    case {1, 4, 5}
        U_fixeddofs = DBoudary.CtrPtsOrd;
        V_fixeddofs = DBoudary.CtrPtsOrd + CtrPts.Num;
    case {2,3}
        U_fixeddofs = DBoudary.CtrPtsOrd1;
        V_fixeddofs = [DBoudary.CtrPtsOrd1; DBoudary.CtrPtsOrd2] + CtrPts.Num;
end
Dofs.Ufixed = U_fixeddofs; Dofs.Vfixed = V_fixeddofs;
Dofs.Free = setdiff(1:Dofs.Num,[Dofs.Ufixed; Dofs.Vfixed]);
U = zeros(Dofs.Num,1);
U(Dofs.Free) = K(Dofs.Free,Dofs.Free)\F(Dofs.Free);
end
%======================================================================================================================%
% Subfunction Solving:                                                                                                 %
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