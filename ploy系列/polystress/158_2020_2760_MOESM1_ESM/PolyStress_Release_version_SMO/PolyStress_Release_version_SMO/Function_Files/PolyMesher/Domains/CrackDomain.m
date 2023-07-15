%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = CrackDomain(Demand,Arg)
  L = 2; % Domain size
  BdBox = [0 L/2 -L/2 L/2];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox,L);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox,L)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  LeftEdgeNodes = find((Node(:,2)<=0)&(abs(Node(:,1))<eps)); % fixed nodes
  BottomNode = find(sqrt(Node(:,1).^2+(Node(:,2)-BdBox(3)).^2)<eps);
  TopNodes = find((Node(:,2)>0.9*L/2-eps)&(abs(Node(:,1)-L/2)<eps));
  FixedNodes = [LeftEdgeNodes];
  Supp = ones(length(FixedNodes),3);
  Supp(:,1) = FixedNodes; 
  Supp(:,3) = 0; Supp(FixedNodes==BottomNode,3) = 1;
  n = length(TopNodes);
  Load = [TopNodes,5*ones(n,1)/n,zeros(n,1)]; % Load array 
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%