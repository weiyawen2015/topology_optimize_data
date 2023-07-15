%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [VolFrac,ElemInd,MatInd,SElemInd,SMatInd,Mat,Color] = ShortColumnConstraints(Seeds)
  NElem = size(Seeds,1);
  NMat = 2;
  NConstr = NMat;
  VolFrac = [0.4; 0.4];
  Mat = [2; 0.3]./2;
  Color = [1 0 0; %red
           0 1 0]; %green
  figure; hold on; box on;
  for i = 1:NMat
    stem(i,Mat(i),'.','Color',Color(i,:),'LineWidth',4);
  end
  xlim([0,NMat+1]); ylim([0,max(Mat)]);
  xlabel('material number'); ylabel('Young''s Modulus, E^0'); 
  set(gca,'XTick',linspace(1,NMat,NMat));
  MatInd = cell(NConstr,1);
  for i = 1:NConstr
      MatInd{i} = i;
  end
  ElemInd = cell(NConstr,1);
  for i = 1:NConstr
    ElemInd{i} = linspace(1,NElem,NElem);
  end
  SElemInd = []; SMatInd = []; %No passive regions