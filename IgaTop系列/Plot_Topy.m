function [X] = Plot_Topy(X, GauPts, CtrPts, DenFied, Pos)
h1 = figure(1); clf; set(h1,'Position',Pos.p1, 'color',[1 1 1]); % the density figure of Control points
plot3(CtrPts.Cordis(1,:),CtrPts.Cordis(2,:),X.CtrPts,'.','color',[0.5 0 0.8]);
axis equal;
h2 = figure(2); clf; set(h2,'Position',Pos.p1,'color',[1 1 1]);  % density figure of Gauss points
plot3(GauPts.PCor(1,:),GauPts.PCor(2,:),X.GauPts,'.','color',[0.5 0 0.8]);
axis equal;
h3 = figure(3); clf; set(h3,'Position',Pos.p2,'color',[1 1 1]);  % DDF
X.DDF = sum(DenFied.N.*X.CtrPts(DenFied.id),2);
X.DDF = reshape(X.DDF,numel(DenFied.U),numel(DenFied.V))';
surf(DenFied.Ux,DenFied.Uy,X.DDF); shading interp; colormap(jet(256)); alpha(0.5);
axis equal; grid off;
h4 = figure(4); clf; set(h4,'Position',Pos.p3,'color',[1 1 1]); 
GauPts_PCor = GauPts.PCor(1:2, X.GauPts>=0.5);
plot(GauPts_PCor(1,:),GauPts_PCor(2,:),'.','color',[0.5 0 0.8]);
axis equal; axis off;
h5 = figure(5); clf; set(h5,'Position',Pos.p4,'color',[1 1 1]);  % structural topology
contourf(DenFied.Ux, DenFied.Uy, X.DDF, [0.5 0.5], 'facecolor', [0.5 0 0.8], 'edgecolor', [1 1 1]); 
axis equal; axis off;
end
%======================================================================================================================%
% Subfunction Plot_Topy:                                                                                               %
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