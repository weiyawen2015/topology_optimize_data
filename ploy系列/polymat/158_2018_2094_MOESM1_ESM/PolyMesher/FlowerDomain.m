%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [x] = FlowerDomain(Demand,Arg) 
  r1 = 1; 
  r2 = 0.125;
  angle = 1.256637061435917;
  LoadPoints = [r1, 0;
                r1*cos(angle), r1*sin(angle);
                r1*cos(2*angle), r1*sin(2*angle);
                r1*cos(3*angle), r1*sin(3*angle);
                r1*cos(4*angle), r1*sin(4*angle)];
  LoadValues = [-sin(0), cos(0);
          -sin(angle), cos(angle);
          -sin(2*angle), cos(2*angle);
          -sin(3*angle), cos(3*angle);
          -sin(4*angle), cos(4*angle)];
  BdBox = [-1 1 -1 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox,r1,r2);
    case('BC');    x = BndryCnds(Arg{:},BdBox,r2,LoadPoints,LoadValues);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox,LoadPoints);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox,r1,r2)
  dOuterCircle = dCircle(P,0,0,r1);
  dInnerCircle = dCircle(P,0,0,r2);
  Dist = dDiff(dOuterCircle,dInnerCircle);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox,r2,LoadPoints,LoadValues)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  for i = 1:size(LoadPoints,1)
      LoadedNodes(i,1) = find(abs(Node(:,1) - LoadPoints(i,1)) < eps & ...
          abs(Node(:,2) - LoadPoints(i,2)) < eps);
  end
  Load = zeros(size(LoadedNodes,1),3);
  Load(:,1) = LoadedNodes; Load(:,2:3) = LoadValues;
  InnerCircleNodes = ...
      find(abs(sqrt(Node(:,1).^2+ Node(:,2).^2)-r2)<eps); 
  Supp = ones(size(InnerCircleNodes,1),3);
  Supp(:,1) = InnerCircleNodes;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox,LoadPoints)
  PFix = LoadPoints;
%-------------------------------------------------------------------------%