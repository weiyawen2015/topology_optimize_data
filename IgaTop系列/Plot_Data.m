function [DenFied, Pos] = Plot_Data(Num, NURBS)
bdwidth = 5; topbdwidth = 30; scnsize = get(0,'ScreenSize'); 
Pos.p1 = [bdwidth, 3/5*scnsize(4)+bdwidth-50, scnsize(3)/2-2*bdwidth, 2*scnsize(4)/5-(topbdwidth+bdwidth)];
Pos.p2 = [Pos.p1(1)+scnsize(3)/2, Pos.p1(2), Pos.p1(3), Pos.p1(4)];
Pos.p3 = [bdwidth, 1/6*scnsize(4)+bdwidth-100, scnsize(3)/2-2*bdwidth, 2*scnsize(4)/5-(topbdwidth + bdwidth)];
Pos.p4 = [Pos.p1(1)+scnsize(3)/2, Pos.p3(2), Pos.p3(3), Pos.p3(4)];
Uknots = linspace(0,1,10*Num(1)); Vknots = linspace(0,1,10*Num(2));
[N_f, id_f] = nrbbasisfun({Uknots, Vknots}, NURBS);
[PCor_U,PCor_W] = nrbeval(NURBS, {Uknots, Vknots});
PCor_U = PCor_U./PCor_W;
PCor_Ux = reshape(PCor_U(1,:),numel(Uknots),numel(Vknots))';
PCor_Uy = reshape(PCor_U(2,:),numel(Uknots),numel(Vknots))';
DenFied.N = N_f; DenFied.id = id_f;
DenFied.U = Uknots; DenFied.V = Vknots;
DenFied.Ux = PCor_Ux; DenFied.Uy = PCor_Uy;
end
%======================================================================================================================%
% Subfunction Plot_Data:                                                                                               %
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