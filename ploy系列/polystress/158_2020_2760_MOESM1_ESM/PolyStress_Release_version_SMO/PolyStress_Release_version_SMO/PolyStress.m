%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [z,V,fem] = PolyStress(fem,opt)
%% Initialize variables
Tol = opt.Tol*(opt.zMax-opt.zMin); TolS = opt.TolS;
z = opt.zIni; P = opt.P; Iter = 0; SM = 10;
lambda = opt.lambda0; mu = opt.mu0; 
BFreq = opt.contB(1); B = opt.contB(2); 
Binc = opt.contB(3); Bmax = opt.contB(4);
[E,dEdy,V,dVdy] = opt.MatIntFnc(P*z,B);
Eid = setdiff((1:fem.NElem)',fem.Passive); % Element ID for active elements
z(fem.Passive) = 1; % Set z = 1 for passive elements
%% Plot initial density distribution and von Mises stress map
[fig_img,hV,hS] = InitialPlot(fem,V,0*V); 
%% MMA parameters
Change = 2*opt.Tol; zold1 = z; zold2 = z; 
L = opt.zMin.*ones(fem.NElem,1); U = opt.zMax.*ones(fem.NElem,1);
AsymInc = opt.AsymInc; AsymDecr = opt.AsymDecr;
%% Optimization iterations
tic; 
while (Iter<opt.MaxIter) && (Change>Tol || max(SM)>1+TolS) %AL steps
  Iter = Iter + 1; 
  for j = 1:opt.MMA_Iter %MMA iterations per AL step
    %% Compute cost functionals and sensitivity analysis
    [~,dJdz,f,h,fem] = AL_Function(fem,E,dEdy,V,dVdy,lambda,mu,P); 
    %% Update design variable and analysis parameters  
    [z(Eid),zold1(Eid),zold2(Eid),L(Eid),U(Eid),AsymInc,AsymDecr,Change]...
        = MMA_unconst(dJdz(Eid),z(Eid),zold1(Eid),zold2(Eid),...
                      L(Eid),U(Eid),Iter,AsymInc,AsymDecr,opt);
    [E,dEdy,V,dVdy] = opt.MatIntFnc(P*z,B);
    SM = E.*fem.VM_Stress0/fem.SLim; % Normalized stress measure
    fprintf(['It:%3i_%1i Obj: %1.3f Max_VM: %1.3f |dJ|: %1.3f ',...
        'Ch/Tol: %1.3f\n'],Iter,j,f,max(SM),norm(dJdz),Change/Tol);
    if (Change<=Tol && max(SM)<=1+TolS), break; end
  end
  %% Update lagrange multiplier estimators and penalty parameter
  lambda = lambda + mu*h; 
  mu = min(opt.alpha*mu,opt.mu_max);
  %% Update material interpolation function
  if mod(Iter,BFreq)==0; B = min(B+Binc,Bmax); end
  %% Update density and stress plots
  set(hV,'FaceColor','flat','CData',1-V); drawnow
  set(hS,'FaceColor','flat','CData',SM); drawnow
   %% ×ögif-------------------------------------------
%    global result index_img
%    aa=getframe(fig_img) ;
%    img =  imresize(aa.cdata,1);
%    result{index_img} = img;   
%    index_img = index_img + 1;
end
t = toc;
if t<=60, fprintf('Optimization time: %i seconds \n', round(t))
elseif t<=3600, fprintf('Optimization time: %1.1f minutes \n', t/60) 
elseif t<=86400, fprintf('Optimization time: %1.1f hours \n', t/3600) 
else, fprintf('Optimization time: %1.1f days \n', t/86400) 
end
%% ----------------------------- NORMALIZED AL FUNCTION AND ITS SENSITIVITY
function [J,dJdz,f,h,fem] = AL_Function(fem,E,dEdy,V,dVdy,lambda,mu,P)
[f,dfdV,dfdE,fem] = ObjectiveFnc(fem,E,V); 
[h,Penal,dPenaldV,dPenaldE,fem] = PenalFnc(fem,E,V,lambda,mu);
N = fem.NElem;
J = f + Penal/N; %Normalized AL function
dJdz = P'*(dEdy.*(dfdE+dPenaldE./N)+dVdy.*(dfdV+dPenaldV./N)); %Sensitivity
%% ----------------------------------- OBJECTIVE FUNCTION (VOLUME FRACTION)
function [f,dfdV,dfdE,fem] = ObjectiveFnc(fem,E,V)
f = sum(fem.ElemArea.*V)/sum(fem.ElemArea); %Mass ratio
dfdV = fem.ElemArea./sum(fem.ElemArea);
dfdE = zeros(size(E));
%% ---------------------------------- PENALTY FUNCTION (STRESS CONSTRAINTS)
function [h,Penal,dPenaldV,dPenaldE,fem] = PenalFnc(fem,E,V,lambda,mu)
[U,fem] = NLFEM(fem,E); %Run nonlinear FEM routine
[fem.VM_Stress0,dVM_dU] = von_Mises_Stress(fem,U); %VM_Stress and Sensit.
dhdVM = zeros(size(V)); dPenaldV = dhdVM; dPenaldE = dhdVM;
s = fem.VM_Stress0/fem.SLim-1;
g = E.*(s.^3+s);    %Polynomial vanishing constraint
h = max(g,-lambda./mu);
Penal = sum(lambda.*h + mu/2.*h.^2); %Penalty term of AL function
a1 = h==g; %Find entries of h==g for sensitivity computation
dhdVM(a1) = E(a1).*(3.*s(a1).^2+1).*1/fem.SLim; %Sensit. of h wrt VM_Stress
dPenaldVM = (lambda+mu.*h).*dhdVM; % Sensit. of penalty term wrt VM_Stress
dPenaldE(a1) = (lambda(a1)+mu.*h(a1)).*(s(a1).^3+s(a1));
% Adjoint method for sensitivity computation
Adjoint_Load = -accumarray(fem.eDof(:),dVM_dU(:).*dPenaldVM(fem.DofE));
Adjoint_Vector(fem.s) = fem.L'\(fem.L\Adjoint_Load(fem.FreeDofs(fem.s),:));
dFintdE = sparse(fem.eDof,fem.DofE,fem.f_NL);
dPenaldE = dPenaldE + (Adjoint_Vector*dFintdE(fem.FreeDofs,:))'; 
%% --------------------------------- VON MISES STRESS AND SENSITIVITY WRT U
function [VM_Stress,dVM_dU] = von_Mises_Stress(fem,U)
V = [1 -1/2 0; -1/2 1 0; 0 0 3]; % von Mises matrix
ElemU = U(fem.eDof);    %ELement displacement vectors
ee_elem = fem.B0*ElemU; %Strains at the cenroid of all elements
ee_elem = reshape(ee_elem,3,[]); %Strains at the cenroid of all elements
[Cauchy_S, D0] = material_model(fem.MatModel,fem.MatParam,ee_elem);
D0 = sparse(fem.rowD,fem.colD,reshape(D0,9*fem.NElem,1));
DB = D0*fem.B0;
VM_Stress = max(sqrt(sum(Cauchy_S.*(V*Cauchy_S))),eps)'; % von Mises stress
dVM_dCauchy = (V*Cauchy_S)./repmat(VM_Stress',3,1);
dVM_dU = DB'*dVM_dCauchy(:); % Sensitivity of VM_Stress wrt U
%% -------------------------------------- MMA UPDATE SCHEME (UNCONSTRAINED)
function [zNew,zold1,zold2,L,U,AsymInc,AsymDecr,Change] = ...
          MMA_unconst(dfdz,z,zold1,zold2,L,U,Iter,AsymInc,AsymDecr,opt)  
zMin = opt.zMin; zMax = opt.zMax; 
move = opt.Move*(zMax-zMin); Osc=opt.Osc; AsymInit = opt.AsymInit;
xmin = max(zMin,z-move); xmax = min(zMax,z+move);
% Compute asymptotes L and U:
AsymInc = min(1+Osc,AsymInc); 
AsymDecr = max(1-2*Osc,AsymDecr);
if Iter<=2
  L = z - AsymInit*(xmax-xmin);  U = z + AsymInit*(xmax-xmin);    
else
  sgn = (z-zold1).*(zold1-zold2);
  s = ones(size(z)); s(sgn>0) = AsymInc; s(sgn<0) = AsymDecr;
  L = z - s.*(zold1 - L); U = z + s.*(U - zold1);
end
% Compute bounds alpha and beta
alpha = 0.9*L + 0.1*z;   beta = 0.9*U + 0.1*z;
alpha = max(xmin,alpha); beta = min(xmax,beta);
% Solve unconstrained subproblem
feps = 0.000001; 
p = (U-z).^2.*(max(dfdz,0)+0.001*abs(dfdz)+feps./(U-L)); 
q = (z-L).^2.*(-min(dfdz,0)+0.001*abs(dfdz)+feps./(U-L));
zCnd = (L.*p - U.*q + (U - L).*sqrt(p.*q))./(p - q);
zNew = max(alpha,min(beta,zCnd));
zold2 = zold1; zold1=z;
Change = sum(abs(zNew-z))/length(z);
%% ----------------------------------------------------------- INITIAL PLOT
function [fig_img ,handle1,handle2] = InitialPlot(fem,z01,z02)
ElemNodes = cellfun(@length,fem.Element); %Number of nodes of each element
Faces = NaN(fem.NElem,max(ElemNodes));    %Populate Faces with NaN
for el = 1:fem.NElem; Faces(el,1:ElemNodes(el)) = fem.Element{el}(:); end
fig_img = figure(5);
ax1 = subplot(1,2,1); title('Element Densities');
patch('Faces',Faces,'Vertices',fem.Node,'FaceVertexCData',0.*z01,...
      'FaceColor','flat','EdgeColor','k','linewidth',1.5);
handle1 = patch('Faces',Faces,'Vertices',fem.Node,'FaceVertexCData',...
                1-z01,'FaceColor','flat','EdgeColor','none');
axis equal; axis off; axis tight; colormap(ax1,gray); caxis([0 1]);
hsp1 = get(gca, 'Position'); % Get position of subplot 1
ax2 = subplot(1,2,2); title('Normalized von Mises Stress')
handle2 = patch('Faces',Faces,'Vertices',fem.Node,'FaceVertexCData',...
                z02,'FaceColor','flat','EdgeColor','none');
axis equal; axis off; axis tight; colormap(ax2,'jet'); c = colorbar;
w = get(c,'Position');
hsp2 = get(gca, 'Position'); % Get position of subplot 2
set(ax1, 'Position', [hsp1(1)-w(3), hsp1(2:end)]);
set(ax2, 'Position', [hsp2(1)-2*w(3), hsp2(2),  hsp1(3:4)]);
drawnow;
%-------------------------------------------------------------------------%