%------------------------------- PolyStress ------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyStress: A Matlab implementation%
% for topology optimization with local stress constraints using the       %
% augmented Lagrangian method", Structural and Multidisciplinary          %
% Optimization, DOI 10.1007/s00158-020-02664-7, 2020                      %
%-------------------------------------------------------------------------%
function [sigma, D] = material_model(MatModel,MatParam,e) 
nE = size(e,2); % Number of elements
switch MatModel
  case 'Bilinear'
    Et = MatParam(1); % Young's modulus in tension 
    Ec = MatParam(2); % Young's modulus in compression
    G = MatParam(3); % Shear modulus
    vt=Et/(2*G)-1; vc=Ec/(2*G)-1;
    Lt = G*(Et-2*G)/(3*G-Et); Lc = G*(Ec-2*G)/(3*G-Ec);
    e33_t = -(Lt/(2*G+Lt))*(e(1,:)+e(2,:)); 
    e33_c = -(Lc/(2*G+Lc))*(e(1,:)+e(2,:));
    g = zeros(nE,1);
    for ii=1:nE
    if e(1,ii)+e(2,ii)+e33_t(ii)>0 
      g(ii)=e(1,ii)+e(2,ii)+e33_t(ii);
    else
      g(ii)=e(1,ii)+e(2,ii)+e33_c(ii);
    end
    end
    E0 = zeros(1,nE); v = E0;
    E0(g<0)=Ec; E0(g>=0)=Et;
    v(g<0)=vc; v(g>=0)=vt; 
    D = zeros(3,3,nE);
    E = E0./(1-v.^2);
    D(1,1,:) = E;    D(1,2,:) = E.*v; D(1,3,:) = 0;
    D(2,1,:) = E.*v; D(2,2,:) = E;    D(2,3,:) = 0;
    D(3,1,:) = 0;    D(3,2,:) = 0;    D(3,3,:) = E.*(1-v)/2; % Plane stress
    sigma = zeros(3,nE);
    sigma(1,:) = squeeze(D(1,1,:))'.*e(1,:)+squeeze(D(1,2,:))'.*e(2,:); 
    sigma(2,:) = squeeze(D(2,1,:))'.*e(1,:)+squeeze(D(2,2,:))'.*e(2,:);
    sigma(3,:) = squeeze(D(3,3,:))'.*e(3,:);
  case 'Ogden'
    mu = MatParam(1:3:end); % Obtain mu from MatParam array
    alpha = MatParam(2:3:end); % Obtain alpha from MatParam array
    beta = MatParam(3:3:end); % Obtain beta from MatParam array
    u = 0.5*(e(1,:)+e(2,:))'; 
    v = sqrt(((e(1,:)-e(2,:))/2).^2+(e(3,:)/2).^2)';
    e1 = (u+v); e2 = (u-v); % Principal strains
    Th = (0.5*atan2(e(3,:),(e(1,:)-e(2,:))))'; % Principal direction
    %%%%%% PLANE STRESS CONDITION %%%%%%
    % We want to find e33 s.t.: s3 = 0;
    e3=0.*e1;
    J = (1+e1).*(1+e2).*(1+e3);
    s3 = sum(mu.*((1+e3).^(alpha-1)-(1+e3).^(-1).*J.^(-alpha.*beta)),2);
    nR0 = norm(s3);
    for ii=1:5 % Use Newton's method to find e3
      J = (1+e1).*(1+e2).*(1+e3); % Jacobian
      s3 = sum(mu.*((1+e3).^(alpha-1)-(1+e3).^(-1).*J.^(-alpha.*beta)),2);
      gs3 = sum(mu.*((alpha-1).*(e3+1).^(alpha-2)+...
               1./(J.^(alpha.*beta).*(e3+1).^2))+(e1+1).*(e2+1).*...
               (alpha.*beta.*mu)./(J.^(alpha.*beta+1).*(e3+1)),2);
      e3 = e3-gs3.^(-1).*s3;
      if norm(s3)<1e-8*nR0; break; end
    end
    %%%%%% PLANE STRESS CONDITION %%%%%%
    J = (1+e1).*(1+e2).*(1+e3); % Update Jacobian with converged e3
    e31 = -sum((1+e1).^(-1).*(1+e3).^(-1).*mu.*alpha.*beta.*...
        J.^(-alpha.*beta),2)./sum((1+e3).^(alpha-2).*...
        mu.*(alpha-1)+(1+e3).^(-2).*mu.*(1+alpha.*beta).*...
        J.^(-alpha.*beta),2); % Derivative of e3 wrt e1
    e32 = -sum((1+e2).^(-1).*(1+e3).^(-1).*mu.*alpha.*beta.*...
        J.^(-alpha.*beta),2)./sum((1+e3).^(alpha-2).*...
        mu.*(alpha-1)+(1+e3).^(-2).*mu.*(1+alpha.*beta).*...
        J.^(-alpha.*beta),2); % Derivative of e3 wrt e2
    s1 = sum(mu.*((1+e1).^(alpha-1)-(1+e1).^(-1).*J.^(-alpha.*beta)),2); 
    s2 = sum(mu.*((1+e2).^(alpha-1)-(1+e2).^(-1).*J.^(-alpha.*beta)),2); 
    s11 = sum(mu.*(alpha-1).*(e1+1).^(alpha-2)+...
            (1+e1).^(-2).*mu.*(1+alpha.*beta).*J.^(-alpha.*beta),2); 
    s22 = sum(mu.*(alpha-1).*(e2+1).^(alpha-2)+...
            (1+e2).^(-2).*mu.*(1+alpha.*beta).*J.^(-alpha.*beta),2); 
    s12 = sum((1+e1).^(-1).*(1+e2).^(-1).*mu.*alpha.*beta.*...
	          J.^(-alpha.*beta),2); 
    s13 = sum((1+e1).^(-1).*(1+e3).^(-1).*mu.*alpha.*beta.*...
	          J.^(-alpha.*beta),2); 
    s23 = sum((1+e2).^(-1).*(1+e3).^(-1).*mu.*alpha.*beta.*...
	          J.^(-alpha.*beta),2); 
    s11 = s11+s13.*e31;
    s12 = s12+s13.*e32;
    s22 = s22+s23.*e32;
    c = cos(Th); s = sin(Th); c4 = cos(4*Th); s4 = sin(4*Th);
    sigma = zeros(3,nE);
    sigma(1,:) = c.^2.*s1+s.^2.*s2;
    sigma(2,:) = s.^2.*s1+c.^2.*s2;
    sigma(3,:) = s.*c.*s1-s.*c.*s2;
    D1=(s1-s2)./(4*(e1-e2));
    D = zeros(3,3,nE);
    idx=e1==e2;
    D(1,1,idx) = s11(idx);
    D(1,2,idx) = s12(idx);
    D(1,3,idx) = 0;
    D(2,1,idx) = D(1,2,idx);
    D(2,2,idx) = s22(idx);
    D(2,3,idx) = 0;
    D(3,1,idx) = 0;
    D(3,2,idx) = 0;
    D(3,3,idx) = 0.5*(s11(idx)-s12(idx));
    D(1,1,~idx) = s11(~idx).*c(~idx).^4 + 2*s12(~idx).*c(~idx).^2.*...
        s(~idx).^2 + s22(~idx).*s(~idx).^4 + D1(~idx).*(1-c4(~idx));
    D(1,2,~idx) = c(~idx).^2.*(s12(~idx).*c(~idx).^2 + s22(~idx).*...
                  s(~idx).^2) + s(~idx).^2.*(s11(~idx).*c(~idx).^2 + ...
		          s12(~idx).*s(~idx).^2) + D1(~idx).*(c4(~idx)-1);
    D(1,3,~idx) = c(~idx).*s(~idx).*(s11(~idx).*c(~idx).^2 + s12(~idx).*...
        s(~idx).^2) - c(~idx).*s(~idx).*(s12(~idx).*c(~idx).^2 + ...
        s22(~idx).*s(~idx).^2) + D1(~idx).*(-s4(~idx));
    D(2,1,~idx) = D(1,2,~idx);
    D(2,2,~idx) = s22(~idx).*c(~idx).^4 + 2.*s12(~idx).*c(~idx).^2.*...
        s(~idx).^2 + s11(~idx).*s(~idx).^4 + D1(~idx).*(1-c4(~idx));
    D(2,3,~idx) = c(~idx).*s(~idx).*(s12(~idx).*c(~idx).^2 + ...
        s11(~idx).*s(~idx).^2) - c(~idx).*s(~idx).*(s22(~idx).*...
        c(~idx).^2 + s12(~idx).*s(~idx).^2) + D1(~idx).*(s4(~idx));
    D(3,1,~idx) = D(1,3,~idx);
    D(3,2,~idx) = D(2,3,~idx);
    D(3,3,~idx) = c(~idx).^2.*s(~idx).^2.*(s11(~idx) - 2.*s12(~idx) + ...
        s22(~idx)) + D1(~idx).*(1+c4(~idx));
otherwise
  error(['Choose a type of material model in the library, ', ...
         'e.g., "Bilinear", "Ogden"'])
end
D = squeeze(D);
%-------------------------------------------------------------------------%