%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function d1 = dPolygon(P,N_p)
% N_p = [x1,y1;x2,y2;...;xN,yN] where xi, yi are the coordinates of each
% node of the polygon 
  d1 = dLine(P,N_p(1,1),N_p(1,2),N_p(2,1),N_p(2,2));
  N_p2=[N_p;N_p(1,:)];
  for n=2:size(N_p,1)
      d2=dLine(P,N_p2(n,1),N_p2(n,2),N_p2(n+1,1),N_p2(n+1,2));
      d1=dIntersect(d2,d1);
  end
%-------------------------------------------------------------------------%  