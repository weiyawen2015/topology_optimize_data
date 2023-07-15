%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
if exist('Global_K_compile_mex','file')==3
  return;
else
  disp('Compiling MEX files...');
  cd('./Function_Files/NLFEM_Code_mex');
  addpath(genpath('./'));
  %% Create configuration object of class 'coder.MexCodeConfig'.
  cfg = coder.config('mex');
  cfg.GenerateReport = true;
  cfg.ReportPotentialDifferences = false;
  %% Define argument types for entry-point 'Global_K_compile'.
  ARGS = cell(1,1);
  ARGS{1} = cell(10,1);
  ARGS{1}{1} = coder.typeof(0,[Inf  1],[1 0]);
  ARG = coder.typeof(0,[1 Inf],[0 1]);
  ARGS{1}{2} = coder.typeof({ARG}, [Inf  1],[1 0]);
  ARGS{1}{3} = coder.typeof('X',[1 Inf],[0 1]);
  ARGS{1}{4} = coder.typeof(0,[1 Inf],[0 1]);
  ARGS{1}{5} = coder.typeof(0);
  ARGS{1}{6} = coder.typeof(0,[Inf  1],[1 0]);
  ARG = coder.typeof(0,[Inf  1],[1 0]);
  ARGS{1}{7} = coder.typeof({ARG}, [10  1],[1 0]);
  ARG = coder.typeof(0,[Inf Inf],[1 1]);
  ARGS{1}{8} = coder.typeof({ARG}, [10  1],[1 0]);
  ARGS{1}{9} = coder.typeof(0,[Inf  2],[1 0]);
  ARGS{1}{10} = coder.typeof(0,[Inf  1],[1 0]);
  %% Invoke MATLAB Coder.
  codegen -config cfg Global_K_compile -args ARGS{1}
  cd('.\..\..\');
end
%-------------------------------------------------------------------------%