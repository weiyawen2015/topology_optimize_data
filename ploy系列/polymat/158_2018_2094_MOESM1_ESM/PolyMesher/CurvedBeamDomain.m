%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [x] = CurvedBeamDomain(Demand,Arg)
  BdBox = [-1 1 0.25 1.5];
  r = 1.5;
  phi = acos(1/r);
  alpha = (pi/2-phi)/2;
  th = alpha+phi;
  x1 = r*cos(th);
  y1 = r*sin(th);
  x2 = 1;
  y2 = sqrt(r^2-x2^2);
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox,r);
    case('BC');    x = BndryCnds(Arg{:},BdBox,r,x1,y1,x2,y2);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox,r,x1,y1);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox,r)
  circle = dCircle(P,0,0,r);
  small_circle = dCircle(P,0,0.25,0.5);
  d = [BdBox(1)-P(:,1), P(:,1)-BdBox(2), BdBox(3)-P(:,2)];
  lines = [d,max(d,[],2)];
  Dist = dIntersect(dDiff(circle,small_circle),lines);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox,r,x1,y1,x2,y2)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  OneThirdNode = find(abs(Node(:,1)+x1)<eps & ...
                       abs(Node(:,2)-y1)<eps);
  MidNode = find(abs(Node(:,1)-0)<eps & ...
                       abs(Node(:,2)-r)<eps);
  TwoThirdNode = find(abs(Node(:,1)-x1)<eps & ...
                       abs(Node(:,2)-y1)<eps);
  LeftNode = find(abs(Node(:,1)+x2)<eps & ...
                       abs(Node(:,2)-y2)<eps);
  RightNode = find(abs(Node(:,1)-x2)<eps & ...
                       abs(Node(:,2)-y2)<eps);
  BottomNodes = find(abs(Node(:,2)-0.25)<eps);
  FixedNodes = BottomNodes;
  Supp = zeros(length(FixedNodes),3);
  Supp(:,1)=FixedNodes; Supp(1:end,2)=1; Supp(1:end,3)=1;
  th1 = sin(y1/r);
  th2 = sin(y2/r);
  Load = [OneThirdNode,1*cos(th1),-1*sin(th1);
          MidNode,0,-1;
          TwoThirdNode,-1*cos(th1),-1*sin(th1);
          LeftNode,0.5*cos(th2),-0.5*sin(th2);
          RightNode,-0.5*cos(th2),-0.5*sin(th2)];
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox,r,x1,y1)
  PFix = [-x1 y1;
          x1 y1;
          0 r];
%-------------------------------------------------------------------------%