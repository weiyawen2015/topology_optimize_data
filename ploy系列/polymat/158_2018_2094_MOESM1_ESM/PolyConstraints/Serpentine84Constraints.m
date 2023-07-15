%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [VolFrac,ElemInd,MatInd,SElemInd,SMatInd,Mat,Color] = Serpentine84Constraints(Seeds)
  % Parameters
  r=4; l=3; a=4;
  b=a*l/r; c=a/r*sqrt(r^2-l^2); d=-sqrt(-l^2+r^2);
  BdBox = [0,3*l+2*b,-r-a-d,r+a+d];
  NMat = 1;
  NConstr = 84;
  VolFrac = 0.5*ones(NConstr,1); 
  % Materials
  [MatInd,Mat,Color] = MaterialsPerConstraint(NConstr,NMat);
  % Constrained regions
  Dist = DistFnc(Seeds,NConstr,BdBox,r,l,a,b,c,d);
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
  for i = 1:NMat
    stem(i,Mat(i),'.','Color',Color(i,:),'LineWidth',4);
  end
  xlim([0,NMat+1]); ylim([0,max(Mat)]);
  xlabel('material number'); ylabel('Young''s Modulus, E^0'); 
  set(gca,'XTick',linspace(1,NMat,NMat));
  % Materials per constraint
  MatInd = cell(NConstr,1);
  for i = 1:size(MatInd,1)
    MatInd{i} = 1;
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,NConstr,BdBox,r,l,a,b,c,d)
Dist = cell(NConstr,1);
circles = cell(14,1); lines = cell(15,1);
%Circles centered at O1 = (0,d)
radius = r+a;
for i = 1:9
    circles{i} = dCircle(P,0,d,radius);
    radius = radius - a/8;
end
%Circles centered at O2 = (2*l+b,c-d)
radius = r+a;
for i = 10:14
    circles{i} = dCircle(P,2*l+b,c-d,radius);
    radius = radius - a/4;
end
%Lines coincident with O1
xyi = [0 d; 0 d]; %Point O1
xyj = [0 1; l 0]; %Endpoints of lines 1 and 9
n = (xyj-xyi)./vecnorm([xyj-xyi],2,2); %unit vectors of lines 1 and 9
angle = acos(dot(n(1,:),n(2,:))/(norm(n(1,:))*norm(n(2,:))))/8;
Rot = [cos(-angle) -sin(-angle); sin(-angle) cos(-angle)];
normal = n(1,:)';
for i = 1:9
    lines{i} = dLine(P,xyi(1,1),xyi(1,2),xyi(1,1)+(r+a)*normal(1),xyi(1,2)+(r+a)*normal(2));
    normal = Rot*normal;
end
%Lines coincident with O2
xyi = [2*l+b c-d; 2*l+b c-d]; %Point O2 
xyj = [l+b,c; 3*l+b,c]; %Endpoints of lines 9 and 15
n = (xyj-xyi)./vecnorm([xyj-xyi],2,2); %unit vectors of lines 9 and 15
angle = acos(dot(n(1,:),n(2,:))/(norm(n(1,:))*norm(n(2,:))))/8;
Rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];
normal = Rot*n(1,:)'; 
for i = 10:13
    lines{i} = -dLine(P,xyi(1,1),xyi(1,2),xyi(1,1)+(r+a)*normal(1),xyi(1,2)+(r+a)*normal(2));
    normal = Rot*normal;
end
normal = Rot*normal; 
for i = 14:15
    lines{i} = -dLine(P,xyi(1,1),xyi(1,2),xyi(1,1)+(r+a)*normal(1),xyi(1,2)+(r+a)*normal(2));
    normal = Rot*(Rot*normal);
end
% Distance values for each constraint
k = 1;
for i = 1:8
    for j = 1:8
        dist = dIntersect(dIntersect(dIntersect(circles{j},-circles{j+1}),-lines{i}),lines{i+1});
        Dist{k} = dist(:,end);
        k = k + 1;
    end
end
for i = 9:12
    for j = 10:13
        dist = dIntersect(dIntersect(dIntersect(circles{j},-circles{j+1}),-lines{i}),lines{i+1});
        Dist{k} = dist(:,end);
        k = k + 1;
    end
end
for i = 13:14
    for j = 10:2:12
        dist = dIntersect(dIntersect(dIntersect(circles{j},-circles{j+2}),-lines{i}),lines{i+1});
        Dist{k} = dist(:,end);
        k = k + 1;
    end
end