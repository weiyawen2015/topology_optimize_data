%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = BridgeDomain(Demand,Arg)
  L = 30; r = L/2; yc = -L/4;
  BdBox = [-L/2 L/2 0 .45*L];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox,r,yc);
    case('BC');    x = BndryCnds(Arg{:},BdBox,L);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox,r,yc)
  circle = dCircle(P,0,yc,r);
  rectangle = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  Dist = dDiff(rectangle,circle);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox,L)
  Lt = 5; % Moving load length
  P = 500e6; % Magnitude of the load
  NStep = 100; % Number of time steps
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  FixedNodes = find(abs(Node(:,1)-BdBox(1))<eps|...
                    abs(Node(:,1)-BdBox(2))<eps);
  DeckNodes = find(abs(Node(:,2)-BdBox(4))<eps); nU=length(DeckNodes);
  Supp = ones(length(FixedNodes),3); Supp(:,1) = FixedNodes;
  Load = zeros(nU,2*(NStep+1)+1);
  Load(:,1) = DeckNodes; 
  for ii=1:NStep+1
    x = L*(ii-1)/NStep-L/2;
    Load((Node(DeckNodes,1)>=x-Lt)&(Node(DeckNodes,1)<=x),2*ii+1) = -P/nU;
  end
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%