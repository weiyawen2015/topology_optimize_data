%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [U,fem] = NLFEM(fem,E)
U = zeros(2*fem.NNode,1); if ~isfield(fem,'U'); fem.U=U; end
if (strcmp(fem.MatModel,'Bilinear')==1)&&(fem.MatParam(1)==fem.MatParam(2))
  K = sparse(fem.i,fem.j,E(fem.e).*fem.k0); % Assemble stiffness matrix
  K = (K+K')/2;
  K = K(fem.FreeDofs,fem.FreeDofs); 
  [U(fem.FreeDofs),L,s] = SolveLinSys(K,fem.Fext(fem.FreeDofs));
  fem.L = L; fem.s = s; % Store Cholesky decomposition information
  fem.f_NL = sparse(fem.iK0,fem.jK0,fem.k0)*U(fem.eDof); 
else
  U = fem.U; % Use previously converged U as initial guess
  [K,~,Res,~,fem] = GlobalK(fem,U,E); %Initial stiffness mtrx. & Res vector
  nRes0 = norm(fem.Fext); nRes = nRes0; % Initial norm of force residual
  % NEWTON-RAPHSON ITERATIONS
  Iter = 0; % Initialize Newton-Raphson iteration counter
  while (nRes>fem.TolR*nRes0 && Iter<=fem.MaxIter)
    [Delta_U, L, s] = SolveLinSys(K,Res);
    [K,~,Res,nRes,fem,U] = LineSearch(U,Delta_U,nRes,-K*Res./nRes,fem,E);
    Iter = Iter+1;
  end
  fem.L = L; fem.s = s; % Store Cholesky decomposition information
  fem.U = U; % Store converged displacement vector
end
%% --------------------------------------------- WEAK LINE SEARCH ALGORITHM 
function [K,Fint,Res,nRes,fem,Un] = LineSearch(U_temp,p,phix,gphix,fem,E)
sigma = 1e-4; alpha_min = 1e-5; alpha_max = 1; MaxIt = 4; 
pgphi = p'*gphix;
alpha = alpha_max;
Un = U_temp; 
Un(fem.FreeDofs) = U_temp(fem.FreeDofs) + alpha*p;
[K,Fint,Res,nRes,fem] = GlobalK(fem,Un,E); phixn=nRes;
it = 0; % Start iteration counter
while (phixn>phix+sigma*alpha*pgphi && alpha>alpha_min && it<=MaxIt) 
  mu = -0.5 * pgphi * alpha / (phixn - phix - alpha * pgphi );
  if mu < 0.1
    mu = 0.5; % Do not trust quadratic interpolation from far away
  end
  alpha = mu*alpha; % New value of line search parameter alpha
  Un(fem.FreeDofs) = U_temp(fem.FreeDofs) + alpha*p; % New disp. vector
  [K,Fint,Res,nRes,fem] = GlobalK(fem,Un,E); phixn=nRes;
  it = it+1;
end
%% ---------------------- GLOBAL STIFFNESS MATRIX AND RESIDUAL FORCE VECTOR
function [K,Fint,Res,nRes,fem] = GlobalK(fem,U,E)
if strcmp(fem.MEX,'Yes')==1
  Compile_mex_file;
  [k0, f_NL] = Global_K_compile_mex(fem.ElemNDof,fem.Element,...
               fem.MatModel,fem.MatParam,fem.Thickness,E,...
               fem.W,fem.dNdxi,fem.Node,U);
elseif strcmp(fem.MEX,'No')==1
  [k0, f_NL] = Global_K_compile(fem.ElemNDof,fem.Element,...
               fem.MatModel,fem.MatParam,fem.Thickness,E,...
               fem.W,fem.dNdxi,fem.Node,U);
else
  error('You should use fem.MEX=Yes or fem.MEX=No');
end
fem.f_NL = f_NL; % Store this for sensitivity computations
f0 = E(fem.DofE).*f_NL;
K = sparse(fem.i,fem.j,E(fem.e).*k0); % Assemble stiffness matrix
K = (K+K')/2; % Symmetrize stiffness matrix
K = K(fem.FreeDofs,fem.FreeDofs);
Fint = sparse(fem.eDof,ones(sum(fem.ElemNDof),1),f0,2*fem.NNode,1);
Res = fem.Fext(fem.FreeDofs)-Fint(fem.FreeDofs); nRes = norm(Res);
%% ----------------------- SOLVE LINEAR SYSTEM USING CHOLESKY DECOMPOSITION
function [Delta_U, L_T, s] = SolveLinSys(K,Res)
[L_T,~,s] = chol(K,'lower','vector');
Delta_U(s,:) = L_T'\(L_T\Res(s,:));
%-------------------------------------------------------------------------%