function [Sh, Hs] = Shep_Fun(CtrPts, rmin)
Ctr_NumU = CtrPts.NumU; Ctr_NumV = CtrPts.NumV;
iH = ones(Ctr_NumU*Ctr_NumV*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH)); sH = zeros(size(iH));
k = 0;
for j1 = 1:Ctr_NumV
    for i1 = 1:Ctr_NumU
        e1 = (j1-1)*Ctr_NumU+i1;
        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),Ctr_NumV)
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),Ctr_NumU)
                e2 = (j2-1)*Ctr_NumU+i2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                theta = sqrt((j1-j2)^2+(i1-i2)^2)./rmin/sqrt(2);
                sH(k) = (max(0, (1-theta)).^6).*(35*theta.^2 + 18*theta + 3);
            end
        end
    end
end
Sh = sparse(iH,jH,sH); Hs = sum(Sh,2);
end
%======================================================================================================================%
% Subfunction Shep_Fun:                                                                                                %
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