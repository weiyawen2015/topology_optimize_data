%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
clear all; clc; close all
restoredefaultpath; addpath(genpath('./')); %Use all folders and subfolders
set(0,'defaulttextinterpreter','latex')
%% ------------------------------------------------------------ CREATE Mesh
[Node,Element,Supp,Load] = Mesh_L_bracket(7000); 
NElem = size(Element,1); % Number of elements
%% ---------------------------------------------------- CREATE 'fem' STRUCT
E0 = 70e3; % E0 in MPa
G = E0/2.5; Et = E0; Ec = E0;  % 0<=(Et,Ec)<=3*G; %Material props. (linear)
fem = struct(...
  'NNode',size(Node,1),...      % Number of nodes
  'NElem',size(Element,1),...   % Number of elements
  'Node',Node,...               % [NNode x 2] array of nodes
  'Element',{Element},...       % [NElement x Var] cell array of elements
  'Supp',Supp,...               % Array of supports
  'Load',Load,...               % Array of loads
  'Passive',[],...              % Passive elements  
  'Thickness',1,...             % Element thickness
  'MatModel','Bilinear',...     % Material model ('Bilinear','Polynomial')
  'MatParam',[Et,Ec,G],...      % Material parameters for MatModel
  'SLim',100,...                % Stress limit
  'TolR', 1e-8, ...             % Tolerance for norm of force residual
  'MaxIter', 15, ...            % Max NR iterations per load step
  'MEX', 'No');                 % Tag to use MEX functions in NLFEM routine
%% ---------------------------------------------------- CREATE 'opt' STRUCT
R = 0.05; q = 3; % Filter radius and filter exponent
p = 3.5; eta0 = 0.5;
m = @(y,B)MatIntFnc(y,'SIMP-H1',[p,B,eta0]);
P = PolyFilter(fem,R,q);
zIni = 0.5*ones(size(P,2),1);
opt = struct(...               
  'zMin',0.0,...              % Lower bound for design variables
  'zMax',1.0,...              % Upper bound for design variables
  'zIni',zIni,...             % Initial design variables
  'MatIntFnc',m,...           % Handle to material interpolation fnc.
  'contB',[5,1,1,10],...      % Threshold projection continuation params.  
  'P',P,...                   % Matrix that maps design to element vars.
  'Tol',0.002,...             % Convergence tolerance on design vars.
  'TolS',0.003,...            % Convergence tolerance on stress constraints
  'MaxIter',150,...           % Maximum number of AL steps
  'MMA_Iter',5,...            % Number of MMA iterations per AL step
  'lambda0',zeros(NElem,1),...% Initial Lagrange multiplier estimators
  'mu0',10,...                % Initial penalty factor for AL function
  'mu_max',10000,...          % Maximum penalty factor for AL function
  'alpha',1.1,...             % Penalty factor update parameter  
  'Move',0.15,...             % Allowable move step in MMA update scheme
  'Osc',0.2,...               % Osc parameter in MMA update scheme
  'AsymInit',0.2,...          % Initial asymptote in MMA update shecme
  'AsymInc',1.2,...           % Asymptote increment in MMA update scheme  
  'AsymDecr',0.7...           % Asymptote decrement in MMA update scheme     
   );
global result index_img
 index_img = 1;
%% ------------------------------------------------------- RUN 'PolyStress'
tic
fem = preComputations(fem); % Run preComputations before running PolyStress
[z,V,fem] = PolyStress(fem,opt);
%-------------------------------------------------------------------------%
toc

% save result result


