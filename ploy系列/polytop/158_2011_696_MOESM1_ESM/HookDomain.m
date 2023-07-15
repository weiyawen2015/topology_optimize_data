%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = HookDomain(Demand,Arg)
  BdBox = [-35.0812, 64.8842, -48.1395, 100.6226];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg);
    case('BC');    x = BndryCnds(Arg,BdBox);
    case('BdBox'); x = BdBox;
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P)
  c1 = dCircle(P,59.9713,78.7683,80);
  c2 = dCircle(P,54.8716,76.8672,35);
  c3 = dCircle(P,0,80.6226,20);
  c4 = dCircle(P,0,80.6226,10);
  c5 = dCircle(P,14.8842,1.8605,50);
  c6 = dCircle(P,0,0,19);
  c7 = dCircle(P,-27.0406,0,8.0406);
  l1 = dLine(P,65.4346,76.9983,-19.9904,81.2407);
  l2 = dLine(P,-25.6060,-27.4746,65.4346,76.9983);
  l3 = dLine(P,1,0,0,0);
  d1 = dDiff(dUnion(dIntersect(dDiff(c1,c2),dIntersect(l1,l2)),c3),c4);
  d2 = dUnion(dIntersect(dDiff(c5,c6),l3),c7);
  d3 = dIntersect(dDiff(c5,c6),-l2);
  Dist = dUnion(dUnion(d1,d2),d3);
  %---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  UpperHalfCircleNodes = find(abs(max(sqrt(Node(:,1).^2+(Node(:,2)...
                                -80.6226).^2)-10,-Node(:,2)+80.6226))<eps);
  Supp = ones(size(UpperHalfCircleNodes,1),3);
  Supp(:,1) = UpperHalfCircleNodes;
  LowerHalfCircleNodes = ...
      find(abs(max(sqrt(Node(:,1).^2+Node(:,2).^2)-19,Node(:,2)))<eps);
  Load = -0.1*ones(size(LowerHalfCircleNodes,1),3);
  Load(:,1) = LowerHalfCircleNodes; Load(:,2) = 0;
  x = {Supp,Load};
%-------------------------------------------------------------------------%