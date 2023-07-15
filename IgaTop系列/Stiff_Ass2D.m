function [K] = Stiff_Ass2D(KE, CtrPts, Ele, Dim, Dofs_Num)
II = zeros(Ele.Num*Dim*Ele.CtrPtsNum*Dim*Ele.CtrPtsNum,1);
JJ = II; KX = II;
ntriplets = 0;
for ide = 1:Ele.Num
    Ele_NoCtPt = Ele.CtrPtsCon(ide,:);
    edof = [Ele_NoCtPt,Ele_NoCtPt+CtrPts.Num];
    for krow = 1:numel(edof)
        for kcol = 1:numel(edof)
            ntriplets = ntriplets+1;
            II(ntriplets) = edof(krow);
            JJ(ntriplets) = edof(kcol);
            KX(ntriplets) = KE{ide}(krow,kcol);
        end
    end
end
K = sparse(II,JJ,KX,Dofs_Num,Dofs_Num); K = (K+K')/2;
end
%======================================================================================================================%
% Subfunction Stiff_Ass2D:                                                                                             %
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