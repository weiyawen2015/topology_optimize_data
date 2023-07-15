%----- ------------------------- PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = EyeBarDomain(Demand,Arg)
  BdBox = [0 1.6 -0.4 0.4];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  d2 = dCircle(P,BdBox(4),0,0.15);
  Dist = dDiff(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  R = 0.15; C = BdBox(4);
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  SupportNodes = find((abs(abs(Node(:,1))-BdBox(2))<eps)&...
                      (abs(Node(:,2))<0.15+eps));
  CircleNodes = find((sqrt((Node(:,1)-C).^2+Node(:,2).^2)<R+eps)&...
                     (Node(:,1)-C<eps));
  yN = Node(CircleNodes,2); P = R^2-yN.^2; P = -70*P./sum(P);  
  Supp = ones(length(SupportNodes),3); Supp(:,1) = SupportNodes;
  n = length(CircleNodes);
  Load = zeros(n,3);
  Load(:,1) = CircleNodes; Load(:,2) = P; 
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%