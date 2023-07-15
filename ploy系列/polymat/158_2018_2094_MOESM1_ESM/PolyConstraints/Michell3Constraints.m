%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [VolFrac,ElemInd,MatInd,SElemInd,SMatInd,Mat,Color] = Michell3Constraints(Seeds)
  % Parameters   
  load('Michell3RegionData.mat');
  BdBox = [0 5 -2 2];
  NMat = 15;
  NConstr = 2*NMat;
  VolFrac = 0.45*ones(NConstr,1);
  % Materials
  [MatInd,Mat,Color] = MaterialsPerConstraint(NConstr,NMat);
  % Constrained regions
  Dist = DistFnc(Seeds,NConstr,RegNode,RegElem,BdBox);
  ElemInd = cell(NConstr,1);
  for i = 1:NConstr
    ElemInd{i} = find(Dist{i}<=0);
  end
  % Passive regions
  SElemInd = []; SMatInd = []; 
%--------------------------------------------- GET MATERIALS PER CONSTRAINT
function [MatInd,Mat,Color] = MaterialsPerConstraint(NConstr,NMat)
  SF = 1000; rng(800);
  Mat = linspace(NMat,1,NMat)'./NMat;
  matorder = randperm(NMat);
  c = jet(SF*NMat);
  % Plot materials
  figure; hold on; box on;
  Color = zeros(NMat,3);
  for i = 1:NMat
    index = matorder(i);
    k = floor((SF*NMat)*Mat(index));
    Color(index,:) = c(k,:);
    stem(index,Mat(index),'.','Color',c(k,:),'LineWidth',4);
  end
  xlabel('material number'); ylabel('Young''s Modulus, E^0'); 
  xlim([0,NMat+1]); ylim([0,max(Mat)]);
  set(gca,'XTick',linspace(1,NMat,NMat));
  grid on; drawnow;
  % Materials per constraint
  MatInd = cell(NConstr,1);
  order = [matorder,matorder];
  for i = 1:NConstr
    MatInd{i} = order(i);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,NConstr,RegNode,RegElem,BdBox)
  d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  d2 = dCircle(P,0,0,BdBox(4)/2);
  Domain = dDiff(d1,d2);
  Dist = cell(NConstr,1);
  for i = 1: NConstr
    d = Domain(:,end);
    Node = cell2mat(RegElem(i));
    NEdge = length(Node);
    for j = 1: NEdge-1
      x1 = RegNode(Node(j),1);
      y1 = RegNode(Node(j),2);
      x2 = RegNode(Node(j+1),1);
      y2 = RegNode(Node(j+1),2);
      d = dIntersect(d,dLine(P,x1,y1,x2,y2));
    end
    x1 = RegNode(Node(NEdge),1);
    y1 = RegNode(Node(NEdge),2);
    x2 = RegNode(Node(1),1);
    y2 = RegNode(Node(1),2);
    d = dIntersect(d,dLine(P,x1,y1,x2,y2));
    Dist{i} = d(:,end);
  end