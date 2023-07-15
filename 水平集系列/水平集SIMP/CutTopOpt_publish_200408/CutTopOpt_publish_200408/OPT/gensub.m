% A MMA IMPLEMENTATION BY NIELS AAGE, JUNE 2012
function [xmma,opt,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
    gensub(k,fx,dfdx,gx,dgdx,x,opt)

%% This function sets up and solve an optimization problem provided as:
%       min_x^n f(x)
%       s.t. g_j(x) < 0,   j = 1,m
%       opt.xmin < x_i < xmax, i = 1,n
%
% INPUT: ALL VECTOR ARE COLUMNS and THE GRADIENT MATRIX IS (N,M)
% k:        integer iteration counter
% dfdx:     (nx1) vector with df/dx1, df/dx2,...,df/dxn
% gx:       (mx1) vector with scalar constraints: g_1, g_2,...,g_m
% dgdx:     (nxm) matrix with dg_1/dx, dg_2/dx,...,dg_m/dx
% opt.U:        (nx1) vector with upper assymptotes
% opt.L:        (nx1) vector with lower assymptotes
% x:        (nx1) vector with x_k, Current iterate

%
% 'opt': is options to control the parameters, heuristics and subproblem solver.
%
% opt.m:      integer size of opt problem
% opt.n:      integer size of opt problem
%
% opt.xold1:      (nx1) vector with x_k-1
% opt.xold2:      (nx1) vector with x_k-2
% opt.xmin:     (nx1) vector with lower bounds
% opt.xmax:     (nx1) vector with upper bounds

% opt.a0:       scalar mma subproblem constant: !!! a0=1 in this code !!!
% opt.a:        (mx1) vector with mma subproblem constants
% opt.c:        (mx1) vector with mma subproblem constants
% opt.d:        (mx1) vector with mma subproblem constants
%
% THE BELOW OPTIONS DOES NOT NEED EXPLICIT SETTING.
% opt.asyminit:         Initial value for k<3 of the asym update (default0.2
% opt.asymdecrease:     Reduction factor for asymptotes (default 0.65)
% opt.asymincrease:     Widening factor for asymptotes (default 1.1)
% opt.constrscaling:    0/1: whether or not the constraint in the
%                            subproblem are artificially scaled (default 0)
%
% OUTPUT:
% xmma:       (nx1) vector with updated design variables
% ymma,...,s: elastic variables, Lagrange multipliers and slacks
% opt:        Updated form of the options struct
%
%% Check input
% more checks are welcome...
if opt.m > opt.n
    error('This implementation is ONLY to be used for n > m')
end

% Check the opt struct and set defaults
if MyIsField(opt,'asyminit')==0
    opt.asyminit=0.15;
    %     opt.asyminit=0.5;
end
if MyIsField(opt,'asymdecrease')==0
    opt.asymdecrease=0.65;
    %     opt.asymdecrease=0.7;
end
if MyIsField(opt,'asymincrease')==0
    opt.asymincrease=1.05;
    %     opt.asymincrease=1.2;
end


if MyIsField(opt,'constrscaling')==0
    opt.constrscaling=0;
end

if length(opt.c)<opt.m
    opt.a=ones(opt.m,1)*opt.a;
    opt.c=ones(opt.m,1)*opt.c;
    opt.d=ones(opt.m,1)*opt.d;
end

%% Set up the mma-subproblem
%
% min  f(x) + a0 z + sum_m(c_i y_i + 0.5 d_i y_i^2)
% s.t. g_i(x) - a_i z - y_i <= 0, i=1,m
% xmin < x_j < xmax, j=1,n
% y_i >= 0, i=1,m
% z >= 0
%
% with
% f   = sum_n(pij./(U-x)+qij./(x-L)), i=1
% gi  = sum_n(pij./(U-x)+qij./(x-L)), i=2,m+1
% pij = see MMA paper 2003
% qij = see MMA paper 2003

% Update assymptotes
if k < 3
    opt.L = x - opt.asyminit*(opt.xmax - opt.xmin);
    opt.U = x + opt.asyminit*(opt.xmax - opt.xmin);
else
    gamma = ones(opt.n,1);
    helpvec = (x-opt.xold1).*(opt.xold1-opt.xold2);
    for i = 1:opt.n
        if helpvec(i) < 0
            gamma(i) = opt.asymdecrease;
        elseif helpvec(i) > 0
            gamma(i) = opt.asymincrease;
        end
    end
    opt.L = x - gamma.*(opt.xold1 - opt.L);
    opt.U = x + gamma.*(opt.U - opt.xold1);
    
    % Robustness correction
    xmi = max(1.0e-5,opt.xmax-opt.xmin);
    opt.L=max(opt.L,x-10.0*xmi);
    opt.L=min(opt.L,x-0.01*xmi);
    opt.U=max(opt.U,x+0.01*xmi);
    opt.U=min(opt.U,x+10*xmi);
    
end

% Update subproblem movelimit
alpha = max(opt.xmin,0.9*opt.L+0.1*x);
beta  = min(opt.xmax,0.9*opt.U+0.1*x);

% Compute coefficients pij and qij for the objective
dfdxp = max(0,dfdx);
dfdxm = max(0,-dfdx);

% Set min number for sens info
feps = 1e-6;
p0 = (opt.U-x).^2.*(dfdxp + 0.001*abs(dfdx) + feps*1./(opt.U-opt.L));
q0 = ((x-opt.L).^2).*(dfdxm + 0.001*abs(dfdx) + feps*1./(opt.U-opt.L));

% For the remaining pij and qij (i = 1:m, j = 1:n)
pij = zeros(opt.n,opt.m);
qij = zeros(opt.n,opt.m);
for i = 1:opt.m
    dgdxp = max(0,dgdx(:,i));
    dgdxm = max(0,-dgdx(:,i));
    if opt.constrscaling==0
        pij(:,i) = ((opt.U-x).^2).*(dgdxp + 0.001*abs(dgdx(:,i)) + feps*1./(opt.U-opt.L));
        qij(:,i) = ((x-opt.L).^2).*(dgdxm + 0.001*abs(dgdx(:,i)) + feps*1./(opt.U-opt.L));
    elseif opt.constrscaling==1
        pij(:,i) = ((opt.U-x).^2).*(dgdxp);
        qij(:,i) = ((x-opt.L).^2).*(dgdxm);
    end
end

% Compute the b_i, i=1:m
b = pij'*(1./(opt.U-x)) +  qij'*(1./(x-opt.L)) - gx;

% Compute the objective correction: !!! NOT USED - just for reference !!!
r = p0'*(1./(opt.U-x)) +  q0'*(1./(x-opt.L)) - fx;

% Solve the subproblem using a dual (primal/dual) interior point method
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
    mmasolve(opt.m,opt.n,opt.L,opt.U,alpha,beta,...
    p0,q0,pij,qij,b,opt.a0,opt.a,opt.c,opt.d);
opt.lam=lam;
end


function isFieldResult = MyIsField(inStruct, fieldName)
% inStruct is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
isFieldResult = 0;
f = fieldnames(inStruct(1));
for i=1:length(f)
    if(strcmp(f{i},strtrim(fieldName)))
        isFieldResult = 1;
        return;
    elseif isstruct(inStruct(1).(f{i}))
        isFieldResult = myIsField(inStruct(1).(f{i}), fieldName);
        if isFieldResult
            return;
        end
    end
end
end
