%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andreasen, C.S., Elingaard, M.O. & Aage, N.                    %
% Level set topology and shape optimization by density methods   %
%    using cut elements with length scale control.               %
% Struct Multidisc Optim (2020).                                 %
% https://doi.org/10.1007/s00158-020-02527-1                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
addpath('FEA')
addpath('MISC')
addpath('GEOM')
addpath('OPT')

ploton =1;
holegen=0;

%% INITIAL STRUCTURE + BC & LOAD
[x,mesh,fem,opt] = init_cantilever();

%% Prepare the geometric constraint & smoothing filter
geom = GeometricConstraintInit(mesh,opt);

%% TOPOPT INIT - set the mma parameters
mma = MMAInit(mesh,opt);

%% Set levelset interval
dx = max(mesh.dx,mesh.dy);
s_min_g = -1/2*dx;
s_max_g = 1/2*dx;

%% Initialize
iter = 0;
change = 1;
gvol = 1;
mmaiter=0;

%% Meshgrid to be used for hole seeding
[XE,YE]=meshgrid(linspace(0+mesh.dx,fem.lx-mesh.dx,fem.nelx),linspace(mesh.dy,fem.ly-mesh.dy,fem.nely));
[XN,YN]=meshgrid(linspace(0,fem.lx,fem.nelx+1),linspace(0,fem.ly,fem.nely+1));

%% Save the desired target volfrac
volfrac_i = opt.volfrac;

% Ensure slow decrease in volume fraction if holegeneration is on
if holegen
    opt.volfrac = 0.95;
end
%% Initialize
mesh1 = mesh; mesh2 = mesh; mesh3 = mesh; % 3 physical realizations
holeinsert=zeros(size(x)); % hole seeding
previdx=[]; % hole seeding
ny=mesh.nely+1;nx=mesh.nelx+1;
xt3old=x;
xt3=x;
change_xt=1;
fval = 1;
history=zeros(opt.max_iter,2);
%% OPTIMIZATION LOOP
while iter < opt.max_iter && (change_xt > 0.01 || fval > 0.0 )
    iter = iter + 1;
    mmaiter = mmaiter +1;
    
    % Filter project
    xf = (geom.Hfilt*x).*geom.Hvec;
    eta1 = 0.5 - opt.deta; % dilated
    eta2 = 0.5 + opt.deta; % eroded
    eta3 = 0.5;   % blueprint
    beta = 12;
    xt3old=xt3;
    xt1 = (tanh(beta*eta1)+tanh(beta*(xf-eta1)))/(tanh(beta*eta1) + tanh(beta*(1-eta1))); % dilated
    xt2 = (tanh(beta*eta2)+tanh(beta*(xf-eta2)))/(tanh(beta*eta2) + tanh(beta*(1-eta2))); % eroded
    xt3 = (tanh(beta*eta3)+tanh(beta*(xf-eta3)))/(tanh(beta*eta3) + tanh(beta*(1-eta3))); % blueprint
    change_xt=norm(xt3-xt3old,'inf');
    
    % Map from x [0; 1] to levelset s [s_min_g s_max_g]
    phi1 = s_min_g + (s_max_g - s_min_g)*xt1; % dilated
    phi2 = s_min_g + (s_max_g - s_min_g)*xt2; % eroded
    phi3 = s_min_g + (s_max_g - s_min_g)*xt3; % Blueprint
    
    % categorize elements and return information to mesh structure
    [mesh1]=categorize_elements(mesh1,phi1);
    [mesh2]=categorize_elements(mesh2,phi2);
    [mesh3]=categorize_elements(mesh3,phi3);
    
    %% FE ANALYSIS + OBJ. FUNC. EVAL. (on eroded design)
    [K2,U2,rho2,mesh2,sedens] = FEA(fem,mesh2);
    
    % objective funcion
    f0val = U2'*K2*U2 ; % objective
    if iter==1
        c_scale=10/f0val;
    end
    
    f0val = f0val*c_scale;
    
    % sensitivity analysis for objective
    [dc] = node_sensitivities(fem,opt,mesh2,phi2,rho2,U2);
    dc = dc*c_scale;
    
    % volume constraint & sensitivity information (on dilated design)
    [curr_vol,dv] = area(mesh1,phi1);
        
    % For showing the area of the blueprint
    [curr_vol3] = area(mesh3,phi3);
    
    % Update volume fraction based on blueprint
    if mod(iter,20)==0 && gvol<0 && ~holegen
        tmp = (volfrac_i * (curr_vol / curr_vol3));
        opt.volfrac = max( tmp, opt.volfrac* 0.95); % make sure not to decrease too fast
        fprintf('Volume constraint updated: opt.volfrac %1.2f (tmp was: %1.2f)\n',opt.volfrac,tmp);
    end
    
    if curr_vol3 <= volfrac_i && holegen == 1
        holegen=0;
        opt.move=0.02;
        fprintf('Disable holegen and set movelimit to %1.2f\n',opt.move);
    end
    
    % Compute constraint and sens.
    gvol= curr_vol/opt.volfrac-1; % volume constraint - on dilate
    dv = dv/opt.volfrac;
    
    % Chain rule for mapping
    df0dx = dc * (s_max_g - s_min_g);
    dv = dv * (s_max_g - s_min_g);
    
    % CHAIN RULE: filter - project
    dx1 = beta * (0.1e1 - tanh(beta * (xf - eta1)) .^ 2) / (tanh(beta * eta1) + tanh(beta * (0.1e1 - eta1)));
    dx2 = beta * (0.1e1 - tanh(beta * (xf - eta2)) .^ 2) / (tanh(beta * eta2) + tanh(beta * (0.1e1 - eta2)));
    
    df0dx = geom.Hfilt*((df0dx.*dx2)./geom.Hfilt_sum);
    dv= geom.Hfilt*((dv.*dx1)./geom.Hfilt_sum);
    
    % Collect constraints
    fval = gvol';
    dfdx=dv; % differentied volume constraint
    
    % Save objective + constraint history
    history(iter,:)=[f0val/c_scale curr_vol3];
    
    %% DESIGN UPDATE
    % update box constraints according to move limit
    mma.xmin=max(0,x-opt.move);
    mma.xmax=min(1,x+opt.move);
    
    % MMA call
    [x_mma,mma] = gensub(mmaiter,f0val,df0dx,fval,dfdx,x,mma);
    change = max(abs(x_mma-x));
    mma.xold2=mma.xold1;
    mma.xold1=x;
    x=x_mma;
    
    %% PLOT DESIGNS USED FOR EVALUATION
    if ploton
        figure(21);clf
        set(gcf,'color','white');
        subplot(2,1,1);
        hold on
        contourf(XN,YN,reshape(phi3,mesh.nely+1,mesh.nelx+1),[0 0],'color','none','facecolor','k');
        contour(XN,YN,reshape(phi1,mesh.nely+1,mesh.nelx+1),[0 0],'r','LineWidth',2);
        contour(XN,YN,reshape(phi3,mesh.nely+1,mesh.nelx+1),[0 0],'k','LineWidth',2);
        contour(XN,YN,reshape(phi2,mesh.nely+1,mesh.nelx+1),[0 0],'m','LineWidth',1);
        axis equal tight
        plot(XN(previdx),YN(previdx),'ro');
        subplot(2,1,2);
        plot(1:iter,history(1:iter,1));
        yyaxis right
        plot(1:iter,history(1:iter,2));
        drawnow
    end
    
    % Evaluate if a hole should be inserted
    sedens_orig=sedens;
    if holegen == 1 && (iter==1 || (mmaiter > 19 ))
        fprintf('Holenucleation condition: curr_vol3 - volfrac_i = %e\n',curr_vol3 - volfrac_i);
        
        % Truncate strain energy density such that no void or cut elemets
        % will be used to introduce a hole.
        sedens(mesh2.void_ele)=max(sedens(:));
        
        % Map to nodes
        rhoNode = zeros((mesh.nely+1)*(mesh.nelx+1),2);
        temp=reshape(sedens,mesh.nely,mesh.nelx);
        elem=0;
        for ex=1:mesh.nelx
            for ey=1:mesh.nely
                elem=elem+1;
                n1 = ey + (mesh.nely+1)*(ex-1);
                n2 = ey + (mesh.nely+1)*ex;
                ndof=[n1 n2 n2+1 n1+1];
                for i=1:4
                    rhoNode(ndof(i),1) = rhoNode(ndof(i),1) + temp(ey,ex);
                    rhoNode(ndof(i),2) = rhoNode(ndof(i),2) + 1;
                end
            end
        end
        sedensn = rhoNode(:,1)./rhoNode(:,2);
        
        %% Find places to seed holes
        tau=(sedensn(:)-min(sedensn(:)))./( max(sedensn(:))-min(sedensn(:))); % map to 0-1
        tau_limit=5e-4;
        idx=find(tau<tau_limit);
        [~,ii]=sort(tau(idx));
        range=1:min(length(ii),10);
        idx=idx(ii(range)); % Only the first 10 candidates holes are seeded
        val=sedensn(idx);
        if ploton
            figure(20);clf;
            
            subplot(2,1,1);
            patch('Faces',mesh.IX,'Vertices', mesh.XY,'FaceVertexCData',log(abs(sedens)'),'FaceColor','flat','edgecolor','none');
            axis equal tight; colorbar,hold on
            plot(mesh.XY(idx,1),mesh.XY(idx,2),'ro');
            title('Truncated ELEMENT strain energy density');
            
            subplot(2,1,2);
            patch(XN(mesh.IX)',YN(mesh.IX)',log(abs(sedensn(mesh.IX))'),'edgecolor','none');
            axis equal tight,  shading interp; colorbar,hold on
            plot(mesh.XY(idx,1),mesh.XY(idx,2),'ro');
            title('Truncated NODAL strain energy density');
        end
        %% Perform hole insertion
        holeidx=idx;
        previdx=[previdx; holeidx;];
        xx = mma.xold1; % use design used for current evaluation to initialize holes
        
        xf1 = (geom.Hfilt*xx).*geom.Hvec; % Normal filter
        xt1 = (tanh(beta*eta1)+tanh(beta*(xf1-eta1)))/(tanh(beta*eta1) + tanh(beta*(1-eta1))); % dilated
        phi_hole = s_min_g + (s_max_g - s_min_g)*xt1; % dilated
        k=0;
        while ~isempty(holeidx) && any(phi_hole(holeidx) > 0) && k < 100 % Until levelset shows a hole
            k=k+1;
            holes=zeros(size(xx));
            holes(holeidx)=0.1;
            xx=max(xx-geom.Hfilt*holes.*geom.Hvec,0); % Perturb design
            
            xf1 = (geom.Hfilt*xx).*geom.Hvec; % Normal filter
            xt1 = (tanh(beta*eta1)+tanh(beta*(xf1-eta1)))/(tanh(beta*eta1) + tanh(beta*(1-eta1))); % dilated
            % Map from x [0; 1] to levelset s [s_min_g s_max_g]
            phi_hole = s_min_g + (s_max_g - s_min_g)*xt1; % dilated
        end
        
        if ~isempty(idx) % Compute volume of design with holes
            [mesh_hole]=categorize_elements(mesh3,phi_hole);
            [curr_vol_hole,~] = area(mesh_hole,phi_hole);
            fprintf('k: %d ',k);
            fprintf('phi_value: %f ',phi_hole(holeidx));
            fprintf('\n');
            x=xx; % update design
            mmaiter = 0;
            
            % updated volume
            tmp = (volfrac_i * (curr_vol / curr_vol3));
            opt.volfrac = max( tmp, opt.volfrac* 0.95); % make sure not to decrease too fast
        end
        
        for i = 1:length(idx)
            fprintf('A hole has been created at: %1.2f %1.2f, sedensn= %1.2e, tau= %1.3e < %1.3e\n',mesh.XY(idx(i),:),val(i),tau(idx(i)),tau_limit);
        end
        
        
    end  
    
    %% print results
    disp([' It.: ' sprintf('%4i',iter)...
        ' Obj.: ' sprintf('%10.4f',full(f0val)) ...
        ' Vol.Const.: ' sprintf('%6.3f',fval(1)) ...
        ' Vol.Dila: ' sprintf('%6.3f',curr_vol )...
        ' Vol.Blue: ' sprintf('%6.3f',curr_vol3 )...
        ' ch.: ' sprintf('%6.3f',change)...
        ' ch.xt:' sprintf('%6.3f',change_xt)...
        ' cut.elem1: ' sprintf('%4i',numel(mesh1.cut_ele)) ...
        ' cut.elem2: ' sprintf('%4i',numel(mesh2.cut_ele))]);
    
end



