%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [VolFrac,ElemInd,MatInd,SElemInd,SMatInd,Mat,Color] = Serpentine1Constraints(Seeds)
  % Parameters  
  VolFrac = 0.5; 
  NMat = 1;
  NConstr = 1;
  % Materials
  Mat = 1;
  MatInd = cell(NMat,1);
  MatInd{1} = 1;
  Color = [0 0 0]; %black
  % Constrained regions
  ElemInd = cell(NConstr,1);
  ElemInd{1} = linspace(1,size(Seeds,1),size(Seeds,1));
  % Passive regions
  SElemInd = []; SMatInd = []; 
  % Plot Materials
  for i = 1:NMat
    stem(i,Mat(i),'.','Color',Color(i,:),'LineWidth',4);
  end
  xlim([0,NMat+1]); ylim([0,max(Mat)]);
  xlabel('material number'); ylabel('Young''s Modulus, E^0'); 
  set(gca,'XTick',linspace(1,NMat,NMat));