%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [Color,RBM] = ConvertColors(Color)
%This function applies a rigid body motion to the RGB cube such that white
%  is at [0 0 0] and black is at [1 1 1]. That way, we plot the
%  multi-material results on a white background.
%Translation Matrix
T1 = eye(4,4); T1(1:3,4) = [-0.5;-0.5;-0.5];
%Rotation Matrix
thz = pi; thx = pi/2;
Rz = [cos(thz) -sin(thz) 0 0;sin(thz) cos(thz) 0 0;0 0 1 0;0 0 0 1];
Rx = [1 0 0 0;0 cos(thx) -sin(thx) 0;0 sin(thx) cos(thx) 0;0 0 0 1];
%Translation Matrix
T2 = eye(4,4); T2(1:3,4) = [0.5;0.5;0.5];
%Full rigid body motion matrix
RBM = T2*Rx*Rz*T1;
%Converted colors
Color = [Color, ones(size(Color,1),1)];
ColorT = RBM*Color';
Color = ColorT';
%-------------------------------------------------------------------------%