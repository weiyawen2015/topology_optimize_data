%-------------------------------- PolyTop --------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyTop: A Matlab %
% implementation of a general topology optimization framework using       %
% unstructured polygonal finite element meshes", Struct Multidisc Optim,  %
% DOI 10.1007/s00158-011-0696-x                                           %
%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [E,dEdy,V,dVdy] = MatIntFnc(y,type,param)
eps = 1e-8;  %Ersatz stiffness
switch(type)
  case('SIMP')
    p = param; % Exponent for SIMP
    V = y; % Volume fraction
    E = eps+(1-eps)*V.^p; % Stiffness interpolation
    dVdy = ones(size(y,1),1);
    dEdy = (1-eps)*p.*y.^(p-1);     
  case('SIMP-H') % SIMP with Heaviside projection 
    p = param(1); % Exponent for SIMP
    B = param(2); % Exponent for Heaviside projection
    V = 1-exp(-B.*y)+y.*exp(-B); % Mass ratio based on heaviside projection
    E = eps+(1-eps)*V.^p; % Stiffness interpolation
    dVdy = B.*exp(-B.*y)+exp(-B); % Derivative of V wrt rho
    dEdy = (1-eps)*p.*V.^(p-1).*dVdy;
  case('SIMP-H1') % SIMP with threshold projection 
    p = param(1); % Exponent for SIMP
    B = param(2); % Penalization parameter for threshold projection
    n = param(3); % Threshold density for threshold projection
    V = (tanh(B*n)+tanh(B*(y-n)))./(tanh(B*n)+tanh(B*(1-n)));
    E = eps+(1-eps)*V.^p;
    dVdy = B*(1-tanh(B*(y-n)).^2)./(tanh(B*n)+tanh(B*(1-n)));
    dEdy = (1-eps)*p*V.^(p-1).*dVdy;
  case('RAMP')
    q = param; % Exponent for RAMP
    V = y;
    E = eps+(1-eps)*V./(1+q*(1-V));
    dVdy = ones(size(y));
    dEdy = ((1-eps)*(q+1))./(q-q*y+1).^2;
  case('RAMP-H')
    q = param(1); % Exponent for RAMP
    B = param(2); % Exponent for Heaviside projection
    V = 1-exp(-B*y)+y*exp(-B);
    E = eps+(1-eps)*V./(1+q*(1-V));
    dVdy = B*exp(-B*y)+exp(-B);
    dEdy = ((1-eps)*(q+1))./(q-q*V+1).^2.*dVdy;
end
%-------------------------------------------------------------------------%