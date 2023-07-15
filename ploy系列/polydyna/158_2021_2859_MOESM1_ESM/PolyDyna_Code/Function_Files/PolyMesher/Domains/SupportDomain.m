%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = SupportDomain(Demand,Arg)
  L = 3;
  BdBox = [-L/2 L/2 0 L];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  NT = 5; % Number of turns of rotating load
  NStep = 150; % Number of time steps
  P = 1e6; % Magnitude of applied load
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  BottomEdgeNodes = find(abs(Node(:,2)-BdBox(3))<eps);
  TopCenterNode = find(abs(Node(:,1))<eps & abs(Node(:,2)-BdBox(4))<eps);
  FixedNodes = BottomEdgeNodes;
  Supp = ones(length(FixedNodes),3);  Supp(:,1)=FixedNodes; 
  Th = linspace(0,NT*2*pi,NStep+1);
  Load = zeros(1,2*(NStep+1)+1); Load(1,1) = TopCenterNode;
  Load(1,2:2:end) = P*cos(Th); Load(1,3:2:end) = P*sin(Th);
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%