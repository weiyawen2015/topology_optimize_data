%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = ClampedDomain(Demand,Arg)
  BdBox = [-6 6 0 2];
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
  NStep = 400; % Number of time steps
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  EdgeNodes = find(abs(Node(:,1)-BdBox(1))<eps | ...
                   abs(Node(:,1)-BdBox(2))<eps);
  BottomCenterNode = find(abs(Node(:,1))<eps & abs(Node(:,2))<eps);
  FixedNodes = EdgeNodes;
  Supp = ones(length(FixedNodes),3);  Supp(:,1)=FixedNodes; 
  Th = linspace(0,pi,NStep+1);
  Load = zeros(1,2*(NStep+1)+1); Load(1,1) = BottomCenterNode;
  Load(1,3:2:end) = -1000*cos(Th); % Half-cycle cosine load
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%