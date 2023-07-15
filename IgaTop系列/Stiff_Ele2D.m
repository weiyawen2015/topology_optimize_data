function [KE, dKE, dv_dg] = Stiff_Ele2D(X, penal, Emin, DH, CtrPts, Ele, GauPts, dRu, dRv)
KE = cell(Ele.Num,1);
dKE = cell(Ele.Num,1);
dv_dg = zeros(GauPts.Num,1);
Nen = Ele.CtrPtsNum;
for ide = 1:Ele.Num
    [idv, idu] = find(Ele.Seque == ide);                    % The two idices in two parametric directions for an element
    Ele_Knot_U = Ele.KnotsU(idu,:);                         % The knot span in the first parametric direction for an element
    Ele_Knot_V = Ele.KnotsV(idv,:);                         % The knot span in the second parametric direction for an element
    Ele_NoCtPt = Ele.CtrPtsCon(ide,:);                      % The number of control points in an element
    Ele_CoCtPt = CtrPts.Cordis(1:2,Ele_NoCtPt);             % The coordinates of the control points in an element
    Ke = zeros(2*Nen,2*Nen);
    dKe = cell(Ele.GauPtsNum,1);
    for i = 1:Ele.GauPtsNum
        GptOrder = GauPts.Seque(ide, i);
        dR_dPara = [dRu(GptOrder,:); dRv(GptOrder,:)];
        dPhy_dPara = dR_dPara*Ele_CoCtPt';
        J1 = dPhy_dPara;
        dR_dPhy = inv(J1)*dR_dPara;
        Be(1,1:Nen) = dR_dPhy(1,:); Be(2,Nen+1:2*Nen) = dR_dPhy(2,:);
        Be(3,1:Nen) = dR_dPhy(2,:); Be(3,Nen+1:2*Nen) = dR_dPhy(1,:);
        dPara_dPare(1,1) = (Ele_Knot_U(2)-Ele_Knot_U(1))/2; % the mapping from the parametric space to the parent space
        dPara_dPare(2,2) = (Ele_Knot_V(2)-Ele_Knot_V(1))/2;
        J2 = dPara_dPare;  J = J1*J2;                       % the mapping from the physical space to the parent space;
        weight = GauPts.Weigh(i)*det(J);                    % Weight factor at this point
        Ke = Ke + (Emin+X.GauPts(GptOrder,:).^penal*(1-Emin))*weight*(Be'*DH*Be);
        dKe{i} = (penal*X.GauPts(GptOrder,:).^(penal-1)*(1-Emin))*weight*(Be'*DH*Be);
        dv_dg(GptOrder) = weight;
    end
    KE{ide} = Ke;
    dKE{ide} = dKe;
end
end
%======================================================================================================================%
% Subfunction Stiff_Ele2D:                                                                                             %
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