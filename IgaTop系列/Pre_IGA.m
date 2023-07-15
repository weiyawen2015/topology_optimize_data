function [CtrPts, Ele, GauPts] = Pre_IGA(NURBS)
%% the unique knots in two parametric directions
Knots.U = unique(NURBS.knots{1})';
Knots.V = unique(NURBS.knots{2})';
%% the information of control points including the physical coordinates, numbers, sequence
CtrPts.Cordis = NURBS.coefs(:,:);
CtrPts.Cordis(1,:) = CtrPts.Cordis(1,:)./CtrPts.Cordis(4,:);   % the X Cartesian coordinates of control points;
CtrPts.Cordis(2,:) = CtrPts.Cordis(2,:)./CtrPts.Cordis(4,:);   % the Y Cartesian coordinates of control points;
CtrPts.Cordis(3,:) = CtrPts.Cordis(3,:)./CtrPts.Cordis(4,:);   % the Z Cartesian coordinates of control points;
CtrPts.Num = prod(NURBS.number);                               % the total number of control points or basis functions;
CtrPts.NumU = NURBS.number(1);                                 % the total number of control points or basis functions in the first parametric;
CtrPts.NumV = NURBS.number(2);                                 % the total number of control points or basis functions in the second parametric;
CtrPts.Seque = reshape(1:CtrPts.Num,CtrPts.NumU,CtrPts.NumV)';
%% the information of the elements (knot spans) in the parametric space, including the numbers, sequence
Ele.NumU = numel(unique(NURBS.knots{1}))-1;                    % the number of elements in the first parametric direction
Ele.NumV = numel(unique(NURBS.knots{2}))-1;                    % the number of elements in the second parametric direction
Ele.Num = Ele.NumU*Ele.NumV;                                   % the number of all elements in the structure
Ele.Seque = reshape(1:Ele.Num, Ele.NumU, Ele.NumV)';
Ele.KnotsU = [Knots.U(1:end-1) Knots.U(2:end)];                % the unique knots of the elements in the first parametric direction
Ele.KnotsV = [Knots.V(1:end-1) Knots.V(2:end)];                % the unique knots of the elements in the second parametric direction
Ele.CtrPtsNum = prod(NURBS.order);
Ele.CtrPtsNumU = NURBS.order(1); Ele.CtrPtsNumV = NURBS.order(2);
[~, Ele.CtrPtsCon] = nrbbasisfun({(sum(Ele.KnotsU,2)./2)', (sum(Ele.KnotsV,2)./2)'}, NURBS);
%% the information of the Gauss quadrature points in the parent space
[GauPts.Weigh, GauPts.QuaPts] = Guadrature(3, numel(NURBS.order));
Ele.GauPtsNum = numel(GauPts.Weigh);
GauPts.Num = Ele.Num*Ele.GauPtsNum;
GauPts.Seque = reshape(1:GauPts.Num,Ele.GauPtsNum,Ele.Num)';
GauPts.CorU = zeros(Ele.Num,Ele.GauPtsNum); GauPts.CorV = zeros(Ele.Num,Ele.GauPtsNum);
for ide = 1:Ele.Num
    [idv, idu] = find(Ele.Seque == ide);                       % The two idices in two parametric directions for an element
    Ele_Knot_U = Ele.KnotsU(idu,:);                            % The knot span in the first parametric direction for an element
    Ele_Knot_V = Ele.KnotsV(idv,:);                            % The knot span in the second parametric direction for an element
    for i = 1:Ele.GauPtsNum
        GauPts.CorU(ide,i) = ((Ele_Knot_U(2)-Ele_Knot_U(1)).*GauPts.QuaPts(i,1) + (Ele_Knot_U(2)+Ele_Knot_U(1)))/2;
        GauPts.CorV(ide,i) = ((Ele_Knot_V(2)-Ele_Knot_V(1)).*GauPts.QuaPts(i,2) + (Ele_Knot_V(2)+Ele_Knot_V(1)))/2;
    end
end
end
%======================================================================================================================%
% Subfunction Pre_IGA:                                                                                                 %
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