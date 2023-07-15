%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [VolFrac,ElemInd,MatInd,SElemInd,SMatInd,Mat,Color] = Flower1MatConstraints(Seeds)
  % Parameters 
  angle = 2*pi/5;
  BdBox = [-1 1 -1 1];
  NMat = 1;
  NConstr = 5;
  VolFrac = 0.3*ones(NConstr,1);
  % Materials
  [MatInd,Mat,Color] = MaterialsPerConstraint(NConstr,NMat);
  % Constrained regions
  Dist = DistFnc(Seeds,NConstr,BdBox,angle);
  ElemInd = cell(NConstr,1);
  for i = 1:NConstr
    ElemInd{i} = find(Dist{i}<=0);
  end
  % Passive regions
  SElemInd = []; SMatInd = [];
  %--------------------------------------------- GET MATERIALS PER CONSTRAINT
function [MatInd,Mat,Color] = MaterialsPerConstraint(NConstr,NMat)
  Mat = 1;
  Color = [0 0 0]; %black
  % Plot materials
  figure; hold on; box on;
  for i = 1:NMat
    stem(i,Mat(i),'.','Color',Color(i,:),'LineWidth',4);
  end
  xlim([0,NMat+1]); ylim([0,max(Mat)]);
  xlabel('material number'); ylabel('Young''s Modulus, E^0'); 
  set(gca,'XTick',linspace(1,NMat,NMat));
  % Materials per constraint
  MatInd = cell(NConstr,1);
  for i = 1:NConstr
      MatInd{i} = 1;
  end
 %----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,NConstr,BdBox,angle)
  % Radial lines  
  d1 = dLine(P,0,0,cos(angle/2),sin(angle/2));
  d2 = dLine (P,0,0,cos(3*angle/2),sin(3*angle/2));
  d3 = dLine(P,0,0,cos(5*angle/2),sin(5*angle/2));
  d4 = dLine(P,0,0,cos(7*angle/2),sin(7*angle/2));
  d5 = dLine(P,0,0,cos(9*angle/2),sin(9*angle/2));
  % Constraints
  g = cell(5,1);
  g{1} = dIntersect(d3,-d4);
  g{2} = dIntersect(d4,-d5);
  g{3} = dIntersect(d5,-d1);
  g{4} = dIntersect(d1,-d2);
  g{5} = dIntersect(d2,-d3);
  for i = 1:NConstr
      Region = g{i};
      Dist{i} = Region(:,end);
  end