%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = SerpentineDomain(Demand,Arg)
  r=4; l=3; a=4;
  b=a*l/r; c=a/r*sqrt(r^2-l^2); d=-sqrt(-l^2+r^2);
  BdBox = [0,3*l+2*b,-r-a-d,r+a+d];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,r,l,a,b,c,d);
    case('BC');    x = BndryCnds(Arg,BdBox,r,l,a,b,c,d);
    case('BdBox'); x = BdBox;
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,r,l,a,b,c,d)
  d1 = dCircle(P,0,d,r+a);
  d2 = dCircle(P,0,d,r);
  d3 = dLine(P,0,d,l,0);
  d4 = dLine(P,0,1,0,d);
  Dist1 = dIntersect(dIntersect(dDiff(d1,d2),d3),d4);
  d5 = dCircle(P,2*l+b,c-d,r+a);
  d6 = dCircle(P,2*l+b,c-d,r);
  d7 = dLine(P,2*l+b,c-d,l+b,c);
  d8 = dLine(P,3*l+b,c,2*l+b,c-d);
  Dist2 = dIntersect(dIntersect(dDiff(d5,d6),d7),d8);
  Dist = dUnion(Dist1,Dist2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,BdBox,r,l,a,b,c,d)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  LeftEdgeNodes = find(abs(Node(:,1)-0)<eps);
  FixedNodes = LeftEdgeNodes;
  Supp = ones(length(FixedNodes),3);
  Supp(:,1)=FixedNodes;
  BottomRightEdge = sqrt((Node(:,1)-(3*l+2*b)).^2+(Node(:,2)).^2);
  [foo,BottomRightEdge] = sort(BottomRightEdge);
  Load = [BottomRightEdge(1),0,-1];
  x = {Supp,Load};
%-------------------------------------------------------------------------%