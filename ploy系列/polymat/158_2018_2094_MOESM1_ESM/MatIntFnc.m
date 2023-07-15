%------------------------------ PolyTop ----------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyTop: A Matlab %
% implementation of a general topology optimization framework using       %
% unstructured polygonal finite element meshes", Struct Multidisc Optim,  %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [w,dwdy,V,dVdy] = MatIntFnc(y,type,param)
switch(type)
  case('SIMP')
    penal = param(1);
    w = y.^penal;
    V = y;
    dwdy = penal*y.^(penal-1);
    dVdy = ones(size(y));
  case('SIMP-H')
    penal = param{1};
    beta = param{2}; 
    h = 1-exp(-beta*y)+y*exp(-beta);
    w = h.^penal;
    V = h;
    dhdy = beta*exp(-beta*y)+exp(-beta);
    dwdy = penal*h.^(penal-1).*dhdy;
    dVdy = dhdy;
  case('RAMP')
    q = param(1);
    w = y./(1+q*(1-y));
    V = y;
    dwdy = (q+1)./(q-q*y+1).^2;
    dVdy = ones(size(y));
  case('RAMP-H')
    q = param{1};
    beta = param(2);
    h = 1-exp(-beta*y)+y*exp(-beta);
    w = h./(1+q*(1-h));
    V = h;
    dhdy = beta*exp(-beta*y)+exp(-beta);
    dwdy = (q+1)./(q-q*h+1).^2.*dhdy;
    dVdy = dhdy;
end
%-------------------------------------------------------------------------%