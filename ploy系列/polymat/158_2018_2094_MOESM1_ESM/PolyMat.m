%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [z,V,fem] = PolyMat(fem,opt,Color)
Iter=0; Tol=opt.Tol*(opt.zMax-opt.zMin); Change=2*Tol*ones(1,opt.NConstr); z=opt.zIni; P=opt.P;
[E,dEdy,V,dVdy] = opt.MultiMatIntFnc(P*z);
[fig_img,FigHandle,FigData,Color,RBM] = InitialPlot(fem,V,Color);
while (Iter<opt.MaxIter) && (max(Change)>Tol)  
  Iter = Iter + 1;
  %Compute cost functionals and analysis sensitivities
  [f,dfdE,dfdV,fem] = ObjectiveFnc(fem,E,V);
  [g,dgdE,dgdV,fem] = ConstraintFnc(fem,opt,E,V);
  %Compute design sensitivities
  dfdz = P'*(dEdy.*repmat(dfdE,1,fem.NMat) + dVdy.*dfdV); 
  for c=1:opt.NConstr
    dgdz(:,c,:) = P'*(dEdy.*repmat(dgdE(:,c),1,fem.NMat) + dVdy.*reshape(dgdV(:,c,:),[fem.NElem,fem.NMat])); 
    %Update design variable and analysis parameters
    ElemIndices = cell2mat(opt.ElemInd(c));
    MatIndices = cell2mat(opt.MatInd(c));
    [z(ElemIndices,MatIndices),Change(c)] = UpdateScheme(...
        dfdz(ElemIndices,MatIndices),g(c),...
        dgdz(ElemIndices,c,MatIndices),z(ElemIndices,MatIndices),V(ElemIndices,MatIndices),opt);
  end
  [E,dEdy,V,dVdy] = opt.MultiMatIntFnc(P*z);
  %Output results
  fprintf('It: %i \t Objective: %1.3f \t Max. Constraint: %1.3f \t Change: %1.3f\n',Iter,f,max(g),max(Change));
  rgb=[zeros(fem.NElem,3),ones(size(V,1),1)];
  for i=1:fem.NMat, rgb(:,1:3) = rgb(:,1:3) + V(:,i)*Color(i,1:3); end
  rgbT = RBM*rgb'; rgb = rgbT';
  I = reshape(rgb(:,1:3),fem.NElem,1,3);
  set(FigHandle,'FaceColor','flat','CData',I(FigData,:,:)); drawnow
     %% ×ögif-------------------------------------------
   global result index_img
   aa=getframe(fig_img) ;
   img =  imresize(aa.cdata,1);
   result{index_img} = img;   
   index_img = index_img + 1;
end
%------------------------------------------------------- OBJECTIVE FUNCTION
function [f,dfdE,dfdV,fem] = ObjectiveFnc(fem,E,V)
[U,fem] = FEAnalysis(fem,E);
f = dot(fem.F,U);
temp = cumsum(-U(fem.i).*fem.k.*U(fem.j));
temp = temp(cumsum(fem.ElemNDof.^2));
dfdE = [temp(1);temp(2:end)-temp(1:end-1)];
dfdV = zeros(size(V));
%----------------------------------------------------- CONSTRAINTS FUNCTION
function [g,dgdE,dgdV,fem] = ConstraintFnc(fem,opt,E,V)
if ~isfield(fem,'ElemArea')
  fem.ElemArea = zeros(fem.NElem,1);
  for el=1:fem.NElem
    vx=fem.Node(fem.Element{el},1); vy=fem.Node(fem.Element{el},2);
    fem.ElemArea(el) = 0.5*sum(vx.*vy([2:end 1])-vy.*vx([2:end 1]));
  end
end
g = zeros(opt.NConstr,1);
dgdV = zeros(fem.NElem,opt.NConstr,fem.NMat);
for c=1:opt.NConstr
  ElemIndices = cell2mat(opt.ElemInd(c));
  MatIndices = cell2mat(opt.MatInd(c));
  for m = 1:size(MatIndices,2)
    g(c) = g(c) + sum(fem.ElemArea(ElemIndices).*V(ElemIndices,MatIndices(m)))./sum(fem.ElemArea(ElemIndices));
    dgdV(ElemIndices,c,MatIndices(m)) = fem.ElemArea(ElemIndices)/sum(fem.ElemArea(ElemIndices));
  end
  g(c) = g(c) - opt.VolFrac(c,1);
end
dgdE = zeros(size(E,1),opt.NConstr); 
%----------------ZPR UPDATE----------------------------------------------- 
function [zNew,Change] = UpdateScheme(dfdz,g,dgdz,z0,V0,opt)
nelem = size(dfdz,1); nmat = size(dfdz,2);
dfdz = reshape(dfdz,nelem*nmat,1); dgdz = reshape(dgdz,nelem*nmat,1);
z0 = reshape(z0,nelem*nmat,1); V0 = reshape(V0,nelem*nmat,1);
zMin=opt.zMin; zMax=opt.zMax;  
move=opt.ZPRMove*(zMax-zMin); eta=opt.ZPREta;
l1=0; l2=1e6;  
while l2-l1 > 1e-4
  lmid = 0.5*(l1+l2);
  B = -(dfdz./dgdz)/lmid;
  zCnd = zMin+(V0-zMin).*max(0,B).^eta; %MultiMatIntFnc may cause non-negative sensitivities
  zNew = max(max(min(min(zCnd,z0+move),zMax),z0-move),zMin);
  if (g+dgdz'*(zNew-z0)>0),  l1=lmid;
  else                       l2=lmid;  end
end
Change = max(abs(zNew-z0))/(zMax-zMin);
zNew = reshape(zNew,nelem,nmat,1);
%-------------------------------------------------------------- FE-ANALYSIS
function [U,fem] = FEAnalysis(fem,E)
if ~isfield(fem,'k')
  fem.ElemNDof = 2*cellfun(@length,fem.Element); % # of DOFs per element
  fem.i = zeros(sum(fem.ElemNDof.^2),1); 
  fem.j=fem.i; fem.k=fem.i; fem.e=fem.i;
  index = 0;
  if ~isfield(fem,'ShapeFnc'), fem=TabShapeFnc(fem); end
  for el = 1:fem.NElem  
    if ~fem.Reg || el==1,  Ke=LocalK(fem,fem.Element{el}); end
    NDof = fem.ElemNDof(el);
    eDof = reshape([2*fem.Element{el}-1;2*fem.Element{el}],NDof,1);
    I=repmat(eDof ,1,NDof); J=I';
    fem.i(index+1:index+NDof^2) = I(:);
    fem.j(index+1:index+NDof^2) = J(:); 
    fem.k(index+1:index+NDof^2) = Ke(:);
    fem.e(index+1:index+NDof^2) = el;
    index = index + NDof^2;
  end
  NLoad = size(fem.Load,1);
  fem.F = zeros(2*fem.NNode,1);  %external load vector
  fem.F(2*fem.Load(1:NLoad,1)-1) = fem.Load(1:NLoad,2);  %x-crdnt
  fem.F(2*fem.Load(1:NLoad,1))   = fem.Load(1:NLoad,3);  %y-crdnt
  NSupp = size(fem.Supp,1);
  FixedDofs = [fem.Supp(1:NSupp,2).*(2*fem.Supp(1:NSupp,1)-1);
               fem.Supp(1:NSupp,3).*(2*fem.Supp(1:NSupp,1))];
  FixedDofs = FixedDofs(FixedDofs>0);
  AllDofs   = [1:2*fem.NNode];
  fem.FreeDofs = setdiff(AllDofs,FixedDofs);
end
K = sparse(fem.i,fem.j,E(fem.e).*fem.k);
K = (K+K')/2;
U = zeros(2*fem.NNode,1);
U(fem.FreeDofs,:) = K(fem.FreeDofs,fem.FreeDofs)\fem.F(fem.FreeDofs,:);
%------------------------------------------------- ELEMENT STIFFNESS MATRIX
function [Ke] = LocalK(fem,eNode)
D=fem.E0/(1-fem.Nu0^2)*[1 fem.Nu0 0;fem.Nu0 1 0;0 0 (1-fem.Nu0)/2];
nn=length(eNode); Ke=zeros(2*nn,2*nn); 
W = fem.ShapeFnc{nn}.W;
for q = 1:length(W)  %quadrature loop
  dNdxi = fem.ShapeFnc{nn}.dNdxi(:,:,q);
  J0 = fem.Node(eNode,:)'*dNdxi; 
  dNdx = dNdxi/J0;
  B = zeros(3,2*nn);
  B(1,1:2:2*nn) = dNdx(:,1)'; 
  B(2,2:2:2*nn) = dNdx(:,2)'; 
  B(3,1:2:2*nn) = dNdx(:,2)'; 
  B(3,2:2:2*nn) = dNdx(:,1)';
  Ke = Ke+B'*D*B*W(q)*det(J0); 
end
%------------------------------------------------- TABULATE SHAPE FUNCTIONS
function fem = TabShapeFnc(fem)
ElemNNode = cellfun(@length,fem.Element); % number of nodes per element
fem.ShapeFnc = cell(max(ElemNNode),1);
for nn = min(ElemNNode):max(ElemNNode)
  [W,Q] = PolyQuad(nn);
  fem.ShapeFnc{nn}.W = W;
  fem.ShapeFnc{nn}.N = zeros(nn,1,size(W,1));
  fem.ShapeFnc{nn}.dNdxi = zeros(nn,2,size(W,1));
  for q = 1:size(W,1)
    [N,dNdxi] = PolyShapeFnc(nn,Q(q,:));
    fem.ShapeFnc{nn}.N(:,:,q) = N;
    fem.ShapeFnc{nn}.dNdxi(:,:,q) = dNdxi;
  end
end
%------------------------------------------------ POLYGONAL SHAPE FUNCTIONS
function [N,dNdxi] = PolyShapeFnc(nn,xi)
N=zeros(nn,1); alpha=zeros(nn,1); dNdxi=zeros(nn,2); dalpha=zeros(nn,2);
sum_alpha=0.0; sum_dalpha=zeros(1,2); A=zeros(nn,1); dA=zeros(nn,2);
[p,Tri] = PolyTrnglt(nn,xi);
for i=1:nn
  sctr = Tri(i,:); pT = p(sctr,:);
  A(i) = 1/2*det([pT,ones(3,1)]);
  dA(i,1) = 1/2*(pT(3,2)-pT(2,2));
  dA(i,2) = 1/2*(pT(2,1)-pT(3,1));
end
A=[A(nn,:);A]; dA=[dA(nn,:);dA];
for i=1:nn
  alpha(i) = 1/(A(i)*A(i+1));
  dalpha(i,1) = -alpha(i)*(dA(i,1)/A(i)+dA(i+1,1)/A(i+1));
  dalpha(i,2) = -alpha(i)*(dA(i,2)/A(i)+dA(i+1,2)/A(i+1));
  sum_alpha = sum_alpha + alpha(i);
  sum_dalpha(1:2) = sum_dalpha(1:2)+dalpha(i,1:2);
end
for i=1:nn
  N(i) = alpha(i)/sum_alpha;
  dNdxi(i,1:2) = (dalpha(i,1:2)-N(i)*sum_dalpha(1:2))/sum_alpha;
end
%---------------------------------------------------- POLYGON TRIANGULATION
function [p,Tri] = PolyTrnglt(nn,xi)
p = [cos(2*pi*((1:nn))/nn); sin(2*pi*((1:nn))/nn)]';
p = [p; xi];
Tri = zeros(nn,3); Tri(1:nn,1)=nn+1;
Tri(1:nn,2)=1:nn; Tri(1:nn,3)=2:nn+1; Tri(nn,3)=1;
%----------------------------------------------------- POLYGONAL QUADRATURE
function [weight,point] = PolyQuad(nn)
[W,Q]= TriQuad;                  %integration pnts & wgts for ref. triangle
[p,Tri] = PolyTrnglt(nn,[0 0]);  %triangulate from origin
point=zeros(nn*length(W),2); weight=zeros(nn*length(W),1);
for k=1:nn
  sctr = Tri(k,:);
  for q=1:length(W)
    [N,dNds] = TriShapeFnc(Q(q,:));  %compute shape functions
    J0 = p(sctr,:)'*dNds;
    l = (k-1)*length(W) + q;
    point(l,:) = N'*p(sctr,:);
    weight(l) = det(J0)*W(q);
  end                                 
end
%---------------------------------------------------- TRIANGULAR QUADRATURE
function [weight,point] = TriQuad
point=[1/6,1/6;2/3,1/6;1/6,2/3]; weight=[1/6,1/6,1/6];   
%----------------------------------------------- TRIANGULAR SHAPE FUNCTIONS
function [N,dNds] = TriShapeFnc(s)
N=[1-s(1)-s(2);s(1);s(2)]; dNds=[-1,-1;1,0;0,1];
%------------------------------------------------------------- INITIAL PLOT
function [fig_img,handle,map,Color,RBM] = InitialPlot(fem,z0,Color)
[Color,RBM] = ConvertColors(Color);
Tri = zeros(length([fem.Element{:}])-2*fem.NElem,3);
map = zeros(size(Tri,1),1); index=0;
for el = 1:fem.NElem
  for enode = 1:length(fem.Element{el})-2
    map(index+1) = el;
    Tri(index+1,:) = fem.Element{el}([1,enode+1,enode+2]);
    index = index + 1;
  end
end
I=[zeros(fem.NElem,3),ones(size(z0,1),1)];
for i=1:fem.NMat
    I(:,1:3) = I(:,1:3) + z0(:,i)*Color(i,1:3);
end
IT = RBM*I';
I = IT';
fig_img = figure(5);
handle = patch('Faces',Tri,'Vertices',fem.Node,'FaceVertexCData',...
               I(map,1:3),'FaceColor','flat','EdgeColor','none');
axis equal; axis off; axis tight; caxis([0 1]);
%-------------------------------------------------------------------------%
%-------------------------- PolyMat - History ----------------------------%
% version: 1.0 (Sept18)
%
% history: Created:    1-Sept-18   Emily Sanders, Anderson Pereira
%          Supervised by:          Glaucio Paulino
%-------------------------------------------------------------------------%