%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function d = dEllipse(P,xc,yc,a,b)
d=a*(sqrt(((P(:,1)-xc)./a).^2+((P(:,2)-yc)./b).^2)-1);
d=[d,d];
%-------------------------------------------------------------------------%