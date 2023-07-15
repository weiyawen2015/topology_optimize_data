%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = L_bracketDomain(Demand,Arg)
  BdBox = [0 1 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),2/5);
  d2 = dRectangle(P,BdBox(1),2/5,BdBox(3),BdBox(4));
  Dist = dUnion(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  TopEdgeNodes = find(abs(Node(:,2)-BdBox(4))<eps);
  RightNodes = find(abs(Node(:,1)-BdBox(2))<eps & ...
                    abs(Node(:,2))>0.85*2/5);  
  FixedNodes = TopEdgeNodes;
  Supp = zeros(length(FixedNodes),3);
  Supp(:,1) = FixedNodes; Supp(1:end,2) = 1; Supp(1:end,3) = 1;
  n = length(RightNodes);
  Load = [RightNodes,zeros(n,1),-2*ones(n,1)/n]; % Load array 
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%