%------------------------------ PolyScript -------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    % 
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
clear; clc; close all
restoredefaultpath; addpath(genpath('./')); %Use all folders and subfolders
set(0,'defaulttextinterpreter','latex')
%% ---------------------------------------------------- CREATE 'fem' STRUCT
Tmax = .05; % Simulation time
[Node,Element,Supp,Load] = Mesh_Cantilever(25088); %close all;
NStep = (size(Load,2)-3)/2;
alpha = 0.05; beta = (1+alpha)^2/4; gamma = (1+2*alpha)/2;%HHT-alpha param.
O = zeros(2*size(Node,1),1);
Obj = 'Compliance';
% Obj = 'Energy';
fem = struct(...
  'NNode',size(Node,1),...     % Number of nodes
  'NElem',size(Element,1),...  % Number of elements
  'Node',Node,...              % [NNode x 2] array of nodes
  'Element',{Element},...      % [NElement x Var] cell array of elements
  'Supp',Supp,...              % Array of supports
  'Load',Load,...              % Array of loads
  'Mass',[],...                % Array of lumped masses
  'u0',O,...                   % Initial displacement vector
  'v0',O,...                   % Initial velocity vector  
  'Thickness',0.01,...         % Element thickness
  'E0',200e9,...               % Young's modulus of solid material
  'Nu0',0.3,...                % Poisson's ratio of solid material  
  'rho',7800,...               % Mass density of solid material (kg/m^3)
  'Ar',[10,1e-5],...           % Rayleigh damping param. C=Ar(1)*M+Ar(2)*K
  'ag',[],...                  % Ground acceleration  
  'Tmax',Tmax,...              % Simulation time
  'NStep',NStep, ...           % Total number of steps  
  'alpha', alpha, ...          % alpha parameter for HHT-alpha method
  'beta', beta, ...            % beta parameter for HHT-alpha method
  'gamma',gamma, ...           % gamma parameter for HHT-alpha method  
  'Obj',Obj,...                % Objective function
  'LL',[],...                  % Vector of DOF index for U_DOF objective
  'Reg',1 ...                  % Tag for regular meshes
   );
%% ---------------------------------------------------- CREATE 'opt' STRUCT
R = 0.2; q = 1;
VolFrac = 0.5;
fem.SElem = []; % Elements in passive solid regions
ElemInd{1} = (1:fem.NElem)'; % Element indices for volume constraint j
P = PolyFilter(fem,R,q,'X');
zIni = ones(size(P,2),1); 
for ii=1:length(VolFrac); zIni(ElemInd{ii})=VolFrac(ii); end % Initial DVs
opt = struct(...               
  'zMin',0,...                 % Lower bound for design variables
  'zMax',1.0,...               % Upper bound for design variables
  'zIni',zIni,...              % Initial design variables
  'MatIntFnc',0,...            % Handle to material interpolation fnc.
  'P',P,...                    % Matrix that maps design to element vars.
  'VolFrac',VolFrac,...        % Arrat if specified volume fraction const.
  'NConstr',size(VolFrac,1),...% Number of volume constraints
  'ElemInd',{ElemInd},...      % Element indices assoc. with each constr.  
  'Tol',0.01,...               % Convergence tolerance on design vars.
  'MaxIter',25,...             % Max. number of optimization iterations
  'Move',0.2,...               % Allowable move step in OC update scheme
  'Eta',0.5 ...                % Exponent used in OC update scheme
   );              
%% ---------------------------------------------------------- RUN 'PolyTop'
fem = PreComputations(fem); % Run preComputations before running PolyStress
figure; B = 1; p_i = 0:1.5:9; tic;
for ii=1:length(p_i)        %Continuation on the RAMP penalty parameter
  if ii==length(p_i); opt.MaxIter = 100; end
  disp(['current p: ', num2str(p_i(ii)), '  current B: ', num2str(B)]);
  opt.MatIntFnc = @(y)MatIntFnc(y,'RAMP-H1',[p_i(ii),B,0.5]);
  [opt.zIni,V,E,fem] = PolyDyna(fem,opt);
  B = min(B+2,10);
end
t0=toc;
if t0<=60, fprintf('Optimization time: %i seconds \n', round(t0))
elseif t0<=3600, fprintf('Optimization time: %1.1f minutes \n', t0/60) 
elseif t0<=86400, fprintf('Optimization time: %1.1f hours \n', t0/3600) 
else, fprintf('Optimization time: %1.1f days \n', t0/86400) 
end
%-------------------------------------------------------------------------%