%------------------------------ PolyTop ----------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyTop: A Matlab %
% implementation of a general topology optimization framework using       %
% unstructured polygonal finite element meshes", Struct Multidisc Optim,  %
% DOI 10.1007/s00158-011-0696-x                                           %
%-------------------------------------------------------------------------%
function [E,dEdy,V,dVdy] = MatIntFnc(y,type,param)
eps = 1e-4;  %Ersatz stiffness
switch(type)
  case('SIMP')
    penal = param;
    E = eps+(1-eps)*y.^penal;
    V = y;
    dEdy = (1-eps)*penal*y.^(penal-1);
    dVdy = ones(size(y,1),1);
  case('SIMP-H')
    penal = param(1);
    beta = param(2); 
    h = 1-exp(-beta*y)+y*exp(-beta);
    E = eps+(1-eps)*h.^penal;
    V = h;
    dhdy = beta*exp(-beta*y)+exp(-beta);
    dEdy = (1-eps)*penal*h.^(penal-1).*dhdy;
    dVdy = dhdy;
  case('RAMP')
    q = param;
    E = eps+(1-eps)*y./(1+q*(1-y));
    V = y;
    dEdy = ((1-eps)*(q+1))./(q-q*y+1).^2;
    dVdy = ones(size(y));
  case('RAMP-H')
    q = param(1);
    beta = param(2);
    h = 1-exp(-beta*y)+y*exp(-beta);
    E = eps+(1-eps)*h./(1+q*(1-h));
    V = h;
    dhdy = beta*exp(-beta*y)+exp(-beta);
    dEdy = ((1-eps)*(q+1))./(q-q*h+1).^2.*dhdy;
    dVdy = dhdy;
end
%-------------------------------------------------------------------------%