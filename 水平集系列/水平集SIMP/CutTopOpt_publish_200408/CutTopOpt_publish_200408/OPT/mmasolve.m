% A MMA SUBPROBLEM SOLVER IMPLEMENTATION BY NIELS AAGE, MAY 2012
function [x,y,z,lam,xi,eta,mu,zeta,s] = ...
    mmasolve(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d)

% function [x,y,z,lam,xi,eta,mu,zet,s] = ...
%    mmasolve(m,n,L,U,alpha,beta,pij,qij,a0,a,b,c,d)
%
% Solves the convex optimization problem:
%
% min_x     sum(p0j./(U-x)+q0j./(x-L)) + a0*z + sum(c.*y + 0.5*d.*y.^2)
%
% s.t.      sum(pij./(U-x)+qij./(x-L)) - ai*z - yi <= bi, i = 1,m
%           alphaj <=  xj <= betaj,  j = 1,n
%           yi >= 0, i = 1,m
%           z >= 0.
%
% NOTE: a0 is fixed at 1 and is not a input parameterin this implementation
%% INITIALIZE
% Lagrange multipliers
lam = c./2.0;

% Barrier
mu = ones(m,1);

% primal vars
x=zeros(n,1);
y=zeros(m,1);
z=0;

% Tolerances and norms
tol=1e-9*sqrt(m+n);
epsi=1.0;
nrI=1.0;
nr2=1.0;

% allocate Newton system
s = zeros(2*m,1);

%% Outer loop for the interior point method
while epsi > tol
    
    loop=0;
    while nrI > 0.9*epsi && loop < 100
        loop=loop+1;
        
        % compute xyz based on lambda
        [x,y,z]=xyzoflambda(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,lam);
        
        % Compute gradient and Hessian
        grad = dualgrad(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,x,y,z,lam,epsi);
        Hess = dualhessian(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,x,y,z,lam,mu);
        
        % solve for search direction
        s(1:m,1) = Hess\grad;
        for i=1:m
            s(m+i,1)= -mu(i)+epsi/lam(i)-s(i)*mu(i)/lam(i);
        end
        
        % Compute linesearch and update lambda and mu
        [lam,mu]=duallinesearch(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,lam,mu,s,epsi);
        
        % New x,y,z
        [x,y,z]=xyzoflambda(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,lam);
        
        % Compute residual at current lambda(x,y,z)
        [nrI,nr2] = dualresidual(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,x,y,z,lam,mu,epsi);
        
        %fprintf('loop: %d, err: %f, search: %f, %f, lam: %f, mu: %f \n',loop,nrI, s(1),s(2),lam(1),mu(1))
        
    end
    % Update the barrier parameter
    epsi=epsi*0.1;
    %fprintf('loop: %d, err: %f, search: %f, %f \n',loop,nrI, s(1),s(2))
    
end


% Just to return something...
xi=[];
eta=[];
mu=[];
zeta=[];


end

function [nrI,nr2] = dualresidual(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,x,y,z,lam,mu,epsi)
% compute the norm of the dual kkt residual
res=zeros(2*m,1);

% compute the constraints
gx = pij'*(1./(U-x)) + qij'*(1./(x-L));

for i=1:m
    res(i)  = -b(i)-y(i)-a(i)*z+gx(i)+mu(i);
    res(i+m)= mu(i)*lam(i) - epsi;
end

%fprintf('res: %f, %f \n',res(1),res(2));

nrI=norm(res,'inf');
nr2=norm(res,2);

end

function [lam,mu]=duallinesearch(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,lam,mu,s,epsi)
% determine stepsize
x0 = [lam; mu];

theta = 1.005;
for i=1:2*m
    if theta < -1.01*s(i)/x0(i)
        theta = -1.01*s(i)/x0(i);
    end
end
step = 1.0/theta;

for i=1:m
    lam(i) = lam(i) + step*s(i);
    mu(i)  = mu(i) + step*s(i+m);
end

end

function [x,y,z]=xyzoflambda(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,lam)
% Compute the x,y and z variables based on lambda
%         // 1.a Correct the lambda's that are less than one
%         // 1.b compute y_i
%         // 1.c compute z
x=zeros(n,1);
y=zeros(m,1);

lamiai=0.0;
for i=1:m
    if (lam(i)<0)
        lam(i)=0.0;
    end
    y(i) = max(0.0,lam(i)-c(i)); % Note y=(lam-c)/d - however d is fixed at one !!
    lamiai = lamiai + lam(i)*a(i);
end
z = max(0.0,10*(lamiai-1)); %// SINCE a0 = 1.0

% 2. compute the x_j
pjlam   = p0 + pij*lam;
qjlam   = q0 + qij*lam;
for j=1:n
    x(j) = (sqrt(pjlam(j))*L(j) + sqrt(qjlam(j))*U(j)) / (sqrt(pjlam(j))+sqrt(qjlam(j)));
    if (x(j)<alpha(j))
        x(j)=alpha(j);
    end
    if (x(j)>beta(j))
        x(j)=beta(j);
    end
end

end

function [grad]=dualgrad(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,x,y,z,lam,epsi)
% Gradient of the dual problem
grad = -1.0*( pij'*(1./(U-x)) +  qij'*(1./(x-L)) - b - a*z - y) - epsi./lam;

end

function [Hess]=dualhessian(m,n,L,U,alpha,beta,p0,q0,pij,qij,b,a0,a,c,d,x,y,z,lam,mu)
% method to compute the hessian of the dual problem + heuristic
% pos.def. check

% help vars
Ux=1./(U-x);
xL=1./(x-L);
Ux2 = Ux.^2;
Ux3 = Ux.^3;
xL2 = xL.^2;
xL3 = xL.^3;

pjlam   = p0 + pij*lam;
qjlam   = q0 + qij*lam;

% hessian wrt x
PQ = pij'*spdiags(Ux2,0,n,n) - qij'*spdiags(xL2,0,n,n);
df2 = 2*pjlam.*Ux3 + 2*qjlam.*xL3; % the 2nd derivative of the objective

% Correct for "out-of-bounds" primal vars.
for j=1:n
    df2(j) = -1.0/df2(j);
    xp(j) = (sqrt(pjlam(j))*L(j) + sqrt(qjlam(j))*U(j)) / (sqrt(pjlam(j))+sqrt(qjlam(j)));
    if (xp(j)<alpha(j))
        df2(j)=0.0;
    end
    if (xp(j)>beta(j))
        df2(j)=0.0;
    end
end

% Hessian contribution from x
Hess = PQ*spdiags(df2,0,n,n)*PQ';

% hessian wrt y and z computation
lamiai=0.0;
for i=1:m
    if lam(i)<0
        lam(i)=0.0;
    end
    lamiai=lamiai+lam(i)*a(i);
    if (lam(i)>c(i))
        Hess(i,i) = Hess(i,i) - 1.0;
    end
    Hess(i,i) = Hess(i,i) - mu(i)/lam(i);
end

% add z contribution if z>=0
if lamiai>1.0
    for i=1:m
        for k=1:m
            Hess(i,k)=Hess(i,k)-10.0*a(i)*a(k);
        end
    end
end

% Ensure pos.def.
HessTrace=0.0;
for i=1:m
    HessTrace=HessTrace+Hess(i,i);
end
HessCorr = 1e-4*HessTrace/m;

% set a lower bound on the correction
if -1.0*HessCorr<1.0e-7
    HessCorr=-1.0e-7;
end

% Add to the hessian
for i=1:m
    Hess(i,i)=Hess(i,i)+HessCorr;
end

end
