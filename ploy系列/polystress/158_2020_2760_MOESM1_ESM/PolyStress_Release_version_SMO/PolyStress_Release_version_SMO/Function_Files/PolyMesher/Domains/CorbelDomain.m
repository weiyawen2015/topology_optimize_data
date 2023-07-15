%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = CorbelDomain(Demand,Arg)
  L = 2;
  BdBox = [0 2*L -1.5*L 1.5*L];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox,L);
    case('BC');    x = BndryCnds(Arg{:},BdBox,L);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox,L)
  d1 = dRectangle(P,BdBox(1),L,BdBox(3),BdBox(4));
  d2 = dRectangle(P,L,BdBox(2),-0.5*L,0.5*L);
  Dist = dUnion(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox,L)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  BottomEdgeNodes = find(abs(Node(:,2)-BdBox(3))<eps);
  TopEdgeNodes = find(abs(Node(:,2)-BdBox(4))<eps);
  RightNodes = find(abs(Node(:,1)-BdBox(2))<eps & ...
      abs(Node(:,2))<0.075*L+eps);    
  FixedNodes = [TopEdgeNodes;BottomEdgeNodes];
  Supp = ones(length(FixedNodes),3);
  Supp(:,1)=FixedNodes; 
  n=length(RightNodes);
  Load = [RightNodes,zeros(n,1),-15*ones(n,1)/n]; % Load array 
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%