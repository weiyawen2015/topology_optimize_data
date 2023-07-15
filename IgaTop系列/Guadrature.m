function [quadweight,quadpoint] = Guadrature(quadorder, dim)
quadpoint  = zeros(quadorder^dim ,dim);
quadweight = zeros(quadorder^dim,1);
r1pt=zeros(quadorder,1);
r1wt=zeros(quadorder,1);
r1pt(1) = 0.774596669241483;
r1pt(2) =-0.774596669241483;
r1pt(3) = 0.000000000000000;
r1wt(1) = 0.555555555555556;
r1wt(2) = 0.555555555555556;
r1wt(3) = 0.888888888888889;
n=1;
for i = 1:quadorder
    for j = 1:quadorder
        quadpoint(n,:) = [ r1pt(i), r1pt(j)];
        quadweight(n)  = r1wt(i)*r1wt(j);
        n = n+1;
    end
end
end
%======================================================================================================================%
% Subfunction Guadrature:                                                                                              %
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