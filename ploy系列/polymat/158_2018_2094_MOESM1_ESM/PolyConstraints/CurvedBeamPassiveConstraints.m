%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [VolFrac,ElemInd,MatInd,SElemInd,SMatInd,Mat,Color] = CurvedBeamPassiveConstraints(Seeds)
  % Parameters 
  BdBox = [-1 1 0.75 1.75];
  r = 1.5;
  x = 0.5; y = sqrt(r^2-x^2);
  NMat = 3;
  NConstr = 3;
  VolFrac = 0.3/NConstr*ones(NConstr,1);
  NFixed = 2;
  % Materials
  [MatInd,SMatInd,Mat,Color] = MaterialsPerConstraint(NConstr,NFixed,NMat);
  % Constrained regions
  [Dist,SDist] = DistFnc(Seeds,NConstr,NFixed,BdBox,r);
  ElemInd = cell(NConstr,1);
  for i = 1:NConstr
    ElemInd{i} = find(Dist{i}<=0);
  end
  % Passive regions
  SElemInd = cell(NFixed,1);
  for i = 1:NFixed
    SElemInd{i} = find(SDist{i}<=0);
  end
  %--------------------------------------------- GET MATERIALS PER CONSTRAINT
function [MatInd,SMatInd,Mat,Color] = MaterialsPerConstraint(NConstr,NFixed,NMat)
  Mat = [1;0.8;0.5];
  Color = [0 0 1; %blue
           0 1 0; %green
           0 1 1];%cyan
  % Plot materials
  figure; hold on; box on;
  for i = 1:NMat
    stem(i,Mat(i),'.','Color',Color(i,:),'LineWidth',4);
  end
  xlim([0,NMat+1]); ylim([0,max(Mat)]);
  xlabel('material number'); ylabel('Young''s Modulus, E^0'); 
  set(gca,'XTick',linspace(1,NMat,NMat));
  % Materials per constraint
  MatInd{1} = [1];
  MatInd{2} = [2];
  MatInd{3} = [3];
  % Materials per passive region
  SMatInd{1} = [1];
  SMatInd{2} = [2];
 %----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function [Dist,SDist] = DistFnc(P,NConstr,NFixed,BdBox,r)
  x1 = 0; y1 = 0; r1 = 1.5; 
  x3 = 0; y3 = 0.25; r3 = 0.5; 
  t = 0.1; 
  c2 = dCircle(P,x1,y1,r1-t);
  c4 = dCircle(P,x3,y3,r3+t);
  dist = dIntersect(c2,-c4);
  Dist = cell(NConstr,1);
  for i = 1:3
    Dist{i} = dist(:,end);
  end
  SDist{1} = -c2(:,end); 
  SDist{2} = c4(:,end);