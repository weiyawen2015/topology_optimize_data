%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [x] = ShortColumnDomain(Demand,Arg) 
  BdBox = [0 1 0 1];
%   LoadedNodes = [1 0.25;
%                  1 0.1875;
%                  1 0.125
%                  1 0.0625
%                  1 0];
%   LoadValues = [1 0;
%                 1 0;
%                 1 0;
%                 1 0;
%                 1 0];
  LoadedNodes = [1 0];
  LoadValues = [1 0];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox,LoadedNodes,LoadValues);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox,LoadedNodes);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  Dist = dRectangle(P,0,1,0,0.5);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox,LoadPoints,LoadValues)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  for i = 1:size(LoadPoints,1)
     LoadedNodes(i,1) = find(abs(Node(:,1) - LoadPoints(i,1)) < eps & ...
          abs(Node(:,2) - LoadPoints(i,2)) < eps);
  end
  Load = zeros(size(LoadedNodes,1),3);
  Load(:,1) = LoadedNodes; Load(:,2:3) = LoadValues;
  FixedNodes = find(abs(Node(:,1)-0)<eps); 
  RollerNodes = find(abs(Node(:,2)-0)<eps);
  SuppFixed = ones(size(FixedNodes,1),3);
  SuppFixed(:,1) = FixedNodes;
  SuppRoller = ones(size(RollerNodes,1),3);
  SuppRoller(:,2) = zeros(size(RollerNodes,1),1);
  SuppRoller(:,1) = RollerNodes;
  Supp = [SuppFixed;SuppRoller];
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox,LoadedNodes)
  PFix = LoadedNodes;
%-------------------------------------------------------------------------%