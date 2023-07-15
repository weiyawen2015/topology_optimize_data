%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = BuildingDomain(Demand,Arg)
  B = 30; H = 75; a = 2/3*H; c = H-a;
  BdBox = [-B/2 B/2 0 H];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox,a,B/2,c);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox,a,b,c)
Th = linspace(0,2*pi,100)';
x = b*cos(Th); y = a*sin(Th)+c;
Ellipse = dPolygon(P,[x,y]);
Rectangle = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
Dist = dIntersect(Rectangle,Ellipse);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  NStep = 80; % Number of time steps
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  FixedNodes = find(abs(Node(:,2)-BdBox(3))<eps);
  Supp = ones(length(FixedNodes),3); Supp(:,1) = FixedNodes;
  TopNode = find(sqrt(Node(:,1).^2+(Node(:,2)-BdBox(4)).^2)<eps);
  Load = zeros(1,2*(NStep+1)+1); Load(1,1) = TopNode;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
PFix = [0,H];
%-------------------------------------------------------------------------%