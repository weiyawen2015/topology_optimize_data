%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = PortalDomain(Demand,Arg)
  BdBox = [-6 6 0 6];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  X1 = [0,3.5;6-0.55,0;6,0;6,6;0,6];
  X2 = [0,3.5;0,6;-6,6;-6,0;-6+0.55,0];
  d1 = dPolygon(P,X1);
  d2 = dPolygon(P,X2);
  Dist = dUnion(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  BottomLeftNodes = find((Node(:,1)<=-5.45+eps)&(abs(Node(:,2))<eps));
  BottomRightNodes = find((Node(:,1)>=5.45-eps)&(abs(Node(:,2))<eps));
  TopNodes = find((abs(Node(:,1))<0.5+eps)&(abs(Node(:,2)-6)<eps));
  Supp = ones(length([BottomLeftNodes;BottomRightNodes]),3);  
  Supp(:,1) = [BottomLeftNodes;BottomRightNodes];
  Supp(length(BottomLeftNodes)+1:end,2) = 0;
  n = length(TopNodes);
  Load = zeros(n,3);
  Load(:,1) = TopNodes; Load(:,3) = -300/n;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%