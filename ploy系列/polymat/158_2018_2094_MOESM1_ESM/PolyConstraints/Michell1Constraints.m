%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [VolFrac,ElemInd,MatInd,SElemInd,SMatInd,Mat,Color] = Michell1Constraints(Seeds)
  % Parameters   
  NElem = size(Seeds,1);
  NMat = 15;
  VolFrac = 0.45;
  % Materials
  [MatInd,Mat,Color] = MaterialsPerConstraint(NMat);
  % Constrained regions
  ElemInd{1} = linspace(1,NElem,NElem);
  % Passive regions
  SElemInd = []; SMatInd = [];
%--------------------------------------------- GET MATERIALS PER CONSTRAINT
function [MatInd,Mat,Color] = MaterialsPerConstraint(NMat)
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
  MatInd{1} = matorder;