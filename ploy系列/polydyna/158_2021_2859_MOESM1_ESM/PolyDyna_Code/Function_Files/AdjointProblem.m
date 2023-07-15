%--------------------------- AdjointProblem ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    %
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
function [xi] = AdjointProblem(fem,M,C,K,dfdu)
al = fem.alpha; B = fem.beta;  g = fem.gamma; % HHT-alpha parameters
DT = fem.Tmax/fem.NStep;  % Time increment
fDOF = fem.FreeDofs;
nu = zeros(2*fem.NNode,fem.NStep+1); xi = nu; mu = nu;
M0 = (1-al)*(1-g)*DT.*C+(1-al)*(1/2-B)*DT^2.*K;
C0 = C+(1-al)*DT*K;
mu(:,end) = dfdu(:,end);
Fext = -B*DT^2*mu(:,end)-g*DT*nu(:,end);
xi(fDOF(fem.s),end) = fem.L'\(fem.L\Fext(fDOF(fem.s)));
for It=fem.NStep+1:-1:3 % Loop over time steps
  mu(:,It-1) = dfdu(:,It-1)+K*xi(:,It)+mu(:,It);
  nu(:,It-1) = C0*xi(:,It)+DT*mu(:,It)+nu(:,It);
  Fext = -M0*xi(:,It)-DT^2*(B*mu(:,It-1)+(1/2-B)*mu(:,It))-...
         DT*(g*nu(:,It-1)+(1-g)*nu(:,It));
  xi(fDOF(fem.s),It-1) = fem.L'\(fem.L\Fext(fDOF(fem.s)));
end
Fext = -M0*xi(:,2)-DT^2*((1/2-B)*mu(:,2))-DT*(1-g)*nu(:,2);  
xi(fDOF,1) = M(fDOF,fDOF)\Fext(fDOF);
%-------------------------------------------------------------------------%