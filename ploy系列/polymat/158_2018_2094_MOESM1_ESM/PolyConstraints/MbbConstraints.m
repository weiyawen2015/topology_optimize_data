%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [VolFrac,ElemInd,MatInd,SElemInd,SMatInd,Mat,Color] = MbbConstraints(Seeds)
BdBox = [0 3 0 1];
VolFrac = 0.5;
NConstr = 1;
NMat = 1;
[MatInd,Mat,Color] = MaterialsPerConstraint(NConstr,NMat);
Dist = DistFnc(Seeds,NConstr,BdBox);
ElemInd = cell(NConstr,1);
for i = 1:NConstr
    ElemInd{i} = find(Dist{i}<=0);
end
SElemInd = []; SMatInd = []; %No passive regions
%--------------------------------------------- GET MATERIALS PER CONSTRAINT
function [MatInd,Mat,Color] = MaterialsPerConstraint(NConstr,NMat)
SF = 1000; rng(1000);
lb = 0.05; ub = 1;
Mat = randperm(NMat)'./NMat;
sortedMat = sort(Mat);
c = jet(SF*NMat);
%%%%% Plot Material Properties with colors
figure; hold on; box on;
Color = zeros(NMat,3);
k = 1;
for i = 1:NMat
    index = find(Mat == sortedMat(i));
    k = floor((SF*NMat)*Mat(index)/ub);
    Color(index,:) = c(k,:);
    stem(index,Mat(index),'.','Color',c(k,:),'LineWidth',4);
    %k = floor(k+SF*NMat/(NMat-1)-1);
end
xlabel('material number'); ylabel('Young''s Modulus, E^0'); 
xlim([0 NMat]); ylim([0 ub]); 
grid on; colormap('jet'); colorbar; drawnow;
%%%%%
MatInd = cell(NConstr,1);
MatInd{1} =[linspace(1,NMat,NMat)];
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,NConstr,BdBox)
Dist = cell(NConstr,1);
Domain = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
for i = 1:NConstr
    Dist{i} = Domain(:,end);
end