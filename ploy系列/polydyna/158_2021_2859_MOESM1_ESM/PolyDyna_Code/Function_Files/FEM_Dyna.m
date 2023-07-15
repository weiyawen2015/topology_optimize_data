%------------------------------ FEM_Dyna ---------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    % 
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
function [un,vn,an,M,C,K,fem] = FEM_Dyna(fem,E,V)
al = fem.alpha; B = fem.beta;  g = fem.gamma; % HHT-alpha parameters
DT = fem.Tmax/fem.NStep;  % Time increment
fDOF = fem.FreeDofs;
an = zeros(2*fem.NNode,fem.NStep+1); un = an; vn = an;
un(:,1) = fem.u0; vn(:,1) = fem.v0; % Initial conditions
[M,C,K] = GlobalMCK(fem,E,V); % Compute stiffness and mass matrices
if ~isempty(fem.Mass)
  DOF = [2*fem.Mass(:,1)-1;2*fem.Mass(:,1)];
  M(DOF,DOF)= M(DOF,DOF) + [fem.Mass(:,2);fem.Mass(:,2)];
end
if ~isempty(fem.ag)
  fem.Fa = -M*ones(2*fem.NNode,1).*fem.ag;
  Ft = fem.Fext+fem.Fa;
else
  Ft = fem.Fext;
end
M1 = M+(1-al)*g*DT.*C+(1-al)*B*DT^2.*K;
an(fDOF,1) = M(fDOF,fDOF)\(Ft(fDOF,1)-K(fDOF,fDOF)*fem.u0(fDOF)...
             -C(fDOF,fDOF)*fem.v0(fDOF)); % Initialize accelerat. vector  
for It=1:fem.NStep % Loop over time steps
  Fext = -C*(vn(:,It)+DT*an(:,It)*(1-al)*(1-g))-...
          K*(an(:,It)*(1/2-B)*(1-al)*DT^2+vn(:,It)*(1-al)*DT+un(:,It))+...
          (1-al)*Ft(:,It+1)+al*Ft(:,It);      
  if It==1
    [an(fDOF,It+1),L,s] = SolveLinSys(M1(fDOF,fDOF),Fext(fDOF));
    fem.L = L; fem.s = s; % Store Cholesky decomposition information
  else
    an(fDOF(s),It+1) = L'\(L\Fext(fDOF(s)));
  end
  un(:,It+1) = un(:,It)+DT*vn(:,It)+(1/2-B)*DT^2*an(:,It)+B*DT^2*an(:,It+1);
  vn(:,It+1) = vn(:,It)+(1-g)*DT*an(:,It)+g*DT*an(:,It+1);
end
%% ------------------------------------- GLOBAL STIFFNESS AND MASS MATRICES
function [M,C,K] = GlobalMCK(fem,E,V)
K = sparse(fem.i,fem.j,E(fem.e).*fem.k0); % Assemble stiffness matrix
K = (K+K')/2; 
M = sparse(fem.i,fem.j,V(fem.e).*fem.m0); % Assemble mass matrix
M = (M+M')/2;
C = fem.Ar(1)*M + fem.Ar(2)*K;            % Compute damping matrix
%% ----------------------- SOLVE LINEAR SYSTEM USING CHOLESKY DECOMPOSITION
function [U, L, s] = SolveLinSys(K,F)
[L,~,s] = chol(K,'lower','vector');
U(s,:) =L'\(L\F(s,:));
%-------------------------------------------------------------------------%