%------------------------------ PolyScript -------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
%% ---------------------------------------------------- CREATE 'fem' STRUCT
close all
clear all
clc
addpath(genpath('./PolyMesher'),genpath('./PolyConstraints'));
[Node,Element,Supp,Load,Seeds] = PolyMesher(@CurvedBeamDomain,10000,300);
[VolFrac,ElemInd,MatInd,SElemInd,SMatInd,Mat,Color] = CurvedBeamConstraints(Seeds);
rmpath('./PolyMesher','./PolyConstraints');
fem = struct(...
  'NNode',size(Node,1),...     % Number of nodes
  'NElem',size(Element,1),...  % Number of elements
  'Node',Node,...              % [NNode x 2] array of nodes
  'Element',{Element},...      % [NElement x Var] cell array of elements
  'Supp',Supp,...              % Array of supports
  'Load',Load,...              % Array of loads
  'Nu0',0.3,...                % Poisson's ratio of solid material
  'E0',1,...                   % Young's modulus of solid material
  'Mat',Mat,...                % Material properties
  'NMat',size(Mat,1),...       % Number of candidate materials
  'Reg',0 ...                  % Tag for regular meshes
   );
%% ---------------------------------------------------- CREATE 'opt' STRUCT
R = 0.03;
P = PolyFilter(fem,R);
zIni = zeros(fem.NElem,fem.NMat);
zIni = InitialGuess(VolFrac,ElemInd,MatInd,SElemInd,SMatInd,zIni);
opt = struct(...               
  'zMin',0.0,...               % Lower bound for design variables
  'zMax',1.0,...               % Upper bound for design variables
  'zIni',zIni,...              % Initial design variables
  'P',P,...                    % Matrix that maps design to element vars.
  'VolFrac',VolFrac,...        % Specified volume fraction constraint
  'NConstr',size(VolFrac,1),...% Number of volume constraints
  'ElemInd',{ElemInd},...      % Element indices assoc. with each constr.
  'MatInd',{MatInd},...        % Material indices assoc. with each constr.
  'Tol',0.01,...               % Convergence tolerance on design vars.
  'MaxIter',100,...            % Max. number of optimization iterations
  'ZPRMove',0.2,...            % Allowable move step in ZPR update scheme
  'ZPREta',1/2 ...             % Exponent used in ZPR update scheme
   );              
%% ----------- RUN 'PolyMat'-----------------------------------------------

Param = [1 1.5 2 3 4;
         0 0.3 0.5 1 1]; %continuation on material interpolation parameters
 global result index_img
 index_img = 1;
for i = 1:size(Param,2)
  penal = Param(1,i); gamma = Param(2,i);
  if i == size(Param,2); opt.P = PolyFilter(fem,-1); end
  disp(['current p: ', num2str(penal), ...
        ', current gamma: ', num2str(gamma)]);
  opt.MatIntFnc = @(y)MatIntFnc(y,'SIMP',penal);
  opt.MultiMatIntFnc = @(y)MultiMatIntFnc(y,opt.MatIntFnc,fem.Mat,gamma);
  [opt.zIni,V,fem] = PolyMat(fem,opt,Color);
end

save result result
%% ------------------------------------------------------------------------