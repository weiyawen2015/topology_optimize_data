%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [E,dEdy,V,dVdy,w] = MultiMatIntFnc(y,MatIntFnc,E0,param)
eps = 1e-4;  %Ersatz stiffness
gamma = param;
[w,dwdy,V,dVdy] = MatIntFnc(y);
NElem = size(y,1); NMat = size(y,2);
dEdw = zeros(NElem,NMat);
% Compute E
S = 1-(gamma.*w); 
Prod = ones(NElem,NMat);
for m = 1:NMat-1
  S = [S(:,NMat),S(:,1:NMat-1)];
  Prod = Prod.*S;
end
E = eps + (1-eps).*((w.*Prod)*E0);
% Compute dEdy
for m = 1:NMat
  S2 = 1-(gamma.*w);
  S2(:,m) = ones(NElem,1);
  dProd = ones(NElem,NMat);
  for j = 1:NMat-1
    S2 = [S2(:,NMat),S2(:,1:NMat-1)];
    dProd = dProd.*S2;
  end  
  w_tmp = -gamma.*w;
  w_tmp(:,m) = ones(NElem,1);
  dEdw(:,m) = (w_tmp.*dProd)*E0;
end
dEdw = dEdw - eps.*dEdw;
dEdy = dwdy.*dEdw;
%-------------------------------------------------------------------------%