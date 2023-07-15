%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = AntennaDomain(Demand,Arg)
  BdBox = [0 3 0 2];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  X1 = [0,0;1,0;0.7,1;0,1];
  X2 = [1,0;3,1.75;2.8,2;0.7,1];
  d1 = dPolygon(P,X1);
  d2 = dPolygon(P,X2);
  Dist = dUnion(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  L = 1; H = 1.75; L1 = 2.8; x1 = 3*L; y1 = H; x2 = L1; y2=2*L;
  d = 0.15*sqrt((x2-x1)^2+(y2-y1)^2);
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  LeftNodes = find(abs(Node(:,1)-BdBox(1))<eps);
  EdgeNodes = find(abs(1.25*(Node(:,1)-2.8)+Node(:,2)-2)<eps&...
                   sqrt((Node(:,1)-2.9).^2+(Node(:,2)-1.875).^2)<d+eps);
  Supp = ones(length(LeftNodes),3);  Supp(:,1)=LeftNodes;
  n = length(EdgeNodes);
  Load = zeros(n,3);
  Load(:,1) = EdgeNodes; Load(:,3) = -40/n;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%