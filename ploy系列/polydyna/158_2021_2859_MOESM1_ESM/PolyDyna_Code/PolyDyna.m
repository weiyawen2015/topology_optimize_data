%------------------------------ PolyDyna ---------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    % 
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
function [z,V,E,fem] = PolyDyna(fem,opt)
Iter=0; Tol=opt.Tol*(opt.zMax-opt.zMin);z=opt.zIni; P=opt.P;
Change = 2*opt.Tol; 
[E,dEdy,V,dVdy] = opt.MatIntFnc(P*z);
[FigHandle] = InitialPlot(fem,V);
c = 0.001*ones(fem.NElem,1); dfdz0=0.*c; % Initialize sens. separ. params
while (Iter<opt.MaxIter) && (max(Change)>Tol) %Optimization iterations
  Iter = Iter + 1;
  %Compute cost functionals and sensitivities
  [f,dfdz,fem] = ObjectiveFnc(fem,E,V,dEdy,dVdy,P);
  [g,dgdz,fem] = ConstraintFnc(fem,opt,E,V,opt.VolFrac,dEdy,dVdy,P);
  %Update design variable and analysis parameters
  for j=1:opt.NConstr
    Eid = cell2mat(opt.ElemInd(j));
    [z(Eid),dfdz0(Eid),c(Eid),Change(j)] = UpdateSchemeSS(...
        dfdz(Eid),dfdz0(Eid),c(Eid),g(j),dgdz(Eid,j),z(Eid),opt);
  end
  [E,dEdy,V,dVdy] = opt.MatIntFnc(P*z);
  %Output results
  fprintf('It: %3i \t Obj.: %1.3f\t Max. Const.: %1.3f\t Change: %1.3f\n',...
           Iter,f,max(g),max(Change));
  set(FigHandle,'FaceColor','flat','CData',1-V); drawnow
end
fprintf('Optimized Objective: %1.3e\n',f/fem.SF);
%------------------------------------------------------- OBJECTIVE FUNCTION
function [f,dfdz,fem] = ObjectiveFnc(fem,E,V,dEdy,dVdy,P)
[U,Ud,Udd,M,C,K,fem] = FEM_Dyna(fem,E,V); %Run dynamics code 
dfdV = zeros(size(V));
dfdE = zeros(size(V));
if strcmp(fem.Obj,'Compliance')==1 %Mean compliance
  f = trace((fem.Fext+fem.Fa)'*U); 
  dfdU = fem.Fext+fem.Fa;
  dfdV = cumsum(sum(fem.Fa0.*U(fem.eDof,:),2));
  dfdV = dfdV(cumsum(fem.ElemNDof)); 
  dfdV = [dfdV(1);dfdV(2:end)-dfdV(1:end-1)]; 
elseif strcmp(fem.Obj,'Energy')==1 %Mean strain energy
  f = 0.5*trace(U'*K*U); 
  dfdU = K*U;
  temp = cumsum(U(fem.i,:).*fem.k0.*U(fem.j,:));
  temp = temp(cumsum(fem.ElemNDof.^2),:);
  dfdE = 0.5*sum([temp(1,:);temp(2:end,:)-temp(1:end-1,:)],2);
elseif strcmp(fem.Obj,'U_DOF')==1 %Mean square displacement at DOF
  f = sum((fem.LL'*U).^2);
  dfdU = 2*(fem.LL'*U).*fem.LL;
end
%Sensitivity analysis using adjoint method
Adj_Vec = AdjointProblem(fem,M,C,K,dfdU); %Solve adjoint problem
al = fem.alpha;
m0 = sparse(fem.iK0,fem.jK0,fem.m0);
k0 = sparse(fem.iK0,fem.jK0,fem.k0);
U = [U(:,1)+fem.Ar(2)*Ud(:,1),(1-al)*(U(:,2:end)+fem.Ar(2)*Ud(:,2:end))+...
                              al*(U(:,1:end-1)+fem.Ar(2)*Ud(:,1:end-1))];
Udd = [Udd(:,1)+fem.Ar(1)*Ud(:,1),Udd(:,2:end)+...
       fem.Ar(1)*((1-al)*Ud(:,2:end)+al*Ud(:,1:end-1))];
Fa0 = [fem.Fa0(:,1), (1-al).*fem.Fa0(:,2:end)+al*fem.Fa0(:,1:end-1)];
for ii=1:fem.NStep+1
  dRdV = m0*Udd(fem.eDof,ii)-Fa0(:,ii);
  dRdV = sparse(fem.eDof,fem.DofE,dRdV);
  dRdE = k0*U(fem.eDof,ii); 
  dRdE = sparse(fem.eDof,fem.DofE,dRdE);
  dfdV = dfdV + (Adj_Vec(fem.FreeDofs,ii)'*dRdV(fem.FreeDofs,:))';
  dfdE = dfdE + (Adj_Vec(fem.FreeDofs,ii)'*dRdE(fem.FreeDofs,:))'; 
end
if ~isfield(fem,'SF'); fem.SF=fem.NStep/f; end %Normalization factor
s = fem.SF;
f = s*f./fem.NStep; dfdV = s*dfdV./fem.NStep; dfdE = s*dfdE./fem.NStep;
dfdz = P'*(dEdy.*dfdE + dVdy.*dfdV); % Chain rule for sensitivity analysis
%------------------------------------------------------ CONSTRAINT FUNCTION
function [g,dgdz,fem] = ConstraintFnc(fem,opt,E,V,VolFrac,dEdy,dVdy,P)
g = zeros(opt.NConstr,1);
dgdV = zeros(fem.NElem,opt.NConstr); dgdE = dgdV;
for j=1:opt.NConstr
  Eid = cell2mat(opt.ElemInd(j));
  g(j) = sum(fem.ElemArea(Eid).*V(Eid))/sum(fem.ElemArea(Eid))-VolFrac(j);
  dgdV(Eid,j) = fem.ElemArea(Eid)/sum(fem.ElemArea(Eid));
end
dgdz = P'*(dEdy.*dgdE + dVdy.*dgdV);
%-------------------------------------------- SENSITIVITY SEPARATION UPDATE
function [zNew,dfdz,c,Change] = UpdateSchemeSS(dfdz,dfdz0,c,g,dgdz,z0,opt)  
zMin=opt.zMin; zMax=opt.zMax;
move=opt.Move*(zMax-zMin); eta=opt.Eta;
dfndz=min(-abs(c).*z0.*eta,dfdz); dfpdz=dfdz-dfndz; %Sensitivity separation
l1=0; l2=1e6;  
while l2-l1 > 1e-4
  lmid = 0.5*(l1+l2);
  B = -dfndz./(dfpdz+lmid*dgdz); 
  zCnd = zMin+(z0-zMin).*B.^eta;
  zNew = max(max(min(min(zCnd,z0+move),zMax),z0-move),zMin);
  if (g+dgdz'*(zNew-z0)>0),  l1=lmid;
  else,                      l2=lmid;  end
end
yk = dfdz-dfdz0; sk = zNew-z0;
al=0.75; zNew = al*zNew+(1-al)*z0;           %Damping scheme
Change = max(abs(zNew-z0))/(zMax-zMin);
c=c+(sk.'*yk-c.'*sk.^2)/(sum(sk.^4)).*sk.^2; %Hessian approximation (PSB)
%------------------------------------------------------------- INITIAL PLOT
function [handle] = InitialPlot(fem,z0)
ElemNodes = cellfun(@length,fem.Element); %Number of nodes of each element
Faces = NaN(fem.NElem,max(ElemNodes));    %Populate Faces with NaN
for el = 1:fem.NElem; Faces(el,1:ElemNodes(el)) = fem.Element{el}(:); end
patch('Faces',Faces,'Vertices',fem.Node,'FaceVertexCData',0.*z0,...
      'FaceColor','flat','EdgeColor','k','linewidth',2);
handle = patch('Faces',Faces,'Vertices',fem.Node,'FaceVertexCData',...
                1-z0,'FaceColor','flat','EdgeColor','none');
axis equal; axis off; axis tight; colormap(gray); caxis([0 1]);
%-------------------------------------------------------------------------%