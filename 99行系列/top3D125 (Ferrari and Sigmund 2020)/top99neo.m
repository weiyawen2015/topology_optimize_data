% function top99neo(nelx,nely,volfrac,penal,rmin,ft,ftBC,eta,beta,move,maxit)
function top99neo
% top99neo(300,100,0.5,3,8.75,3,'N',0.5,2,0.2,500 );
nelx = 120;
nely = 40;
volfrac=0.5;
penal=3;
rmin=3.75;
ft=3;
ftBC = 'N';
eta = 0.5;
beta = 2;
move = 0.2;
maxit = 100;

% ---------------------------- PRE. 1) MATERIAL AND CONTINUATION PARAMETERS
E0 = 1;                                                                    % Young modulus of solid
Emin = 1e-9;                                                               % Young modulus of "void"
nu = 0.3;                                                                  % Poisson ratio
penalCnt = { 1,  1, 25, 0.25 };                                            % continuation scheme on penal
betaCnt  = { 1,  1, 25,    2 };                                            % continuation scheme on beta
if ftBC == 'N', bcF = 'symmetric'; else, bcF = 0; end                      % filter BC selector
% ----------------------------------------- PRE. 2) DISCRETIZATION FEATURES
nEl = nelx * nely;                                                         % number of elements
nodeNrs = int32( reshape( 1 : (1 + nelx) * (1 + nely), 1+nely, 1+nelx ) ); % nodes numbers (defined as int32)
cVec = reshape( 2 * nodeNrs( 1 : end - 1, 1 : end - 1 ) + 1, nEl, 1 );
cMat = cVec + int32( [ 0, 1, 2 * nely + [ 2, 3, 0, 1 ], -2, -1 ] );        % connectivity matrix
nDof = ( 1 + nely ) * ( 1 + nelx ) * 2;                                    % total number of DOFs
[ sI, sII ] = deal( [ ] );
for j = 1 : 8
    sI = cat( 2, sI, j : 8 );
    sII = cat( 2, sII, repmat( j, 1, 8 - j + 1 ) );
end
[ iK , jK ] = deal( cMat( :,  sI )', cMat( :, sII )' );
Iar = sort( [ iK( : ), jK( : ) ], 2, 'descend' ); clear iK jK              % reduced assembly indexing
c1 = [12;3;-6;-3;-6;-3;0;3;12;3;0;-3;-6;-3;-6;12;-3;0;-3;-6;3;12;3;...
    -6;3;-6;12;3;-6;-3;12;3;0;12;-3;12];
c2 = [-4;3;-2;9;2;-3;4;-9;-4;-9;4;-3;2;9;-2;-4;-3;4;9;2;3;-4;-9;-2;...
    3;2;-4;3;-2;9;-4;-9;4;-4;-3;-4];
Ke = 1/(1-nu^2)/24*( c1 + nu .* c2 );                                      % lower sym. part of el. matrix
Ke0( tril( ones( 8 ) ) == 1 ) = Ke';
Ke0 = reshape( Ke0, 8, 8 );
Ke0 = Ke0 + Ke0' - diag( diag( Ke0 ) );                                    % recover full elemental matrix
% ----------------------------- PRE. 3) LOADS, SUPPORTS AND PASSIVE DOMAINS
lcDof = 2 * nodeNrs( 1, 1 );                                               % DOFs with applied load
fixed = union( 1 : 2 : 2*( nely + 1 ), 2 * nodeNrs( end, end ) );          % restrained DOFs
[ pasS, pasV ] = deal( [], [] );                                           % UD, passive solid and void el.
F = fsparse( lcDof', 1, -1, [ nDof, 1 ] );                                 % define load vector
free = setdiff( 1 : nDof, fixed );                                         % set of free DOFs
act = setdiff( (1 : nEl )', union( pasS, pasV ) );                         % set of active d.v.
% --------------------------------------- PRE. 4) DEFINE IMPLICIT FUNCTIONS
prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v(:)-eta)))./...
    (tanh(beta*eta)+tanh(beta*(1-eta)));                                   % projection
deta = @(v,eta,beta) - beta * csch( beta ) .* sech( beta * ( v( : ) - eta ) ).^2 .* ...
    sinh( v( : ) * beta ) .* sinh( ( 1 - v( : ) ) * beta );                % projection eta-derivative
dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));% proj. x-derivative
cnt = @(v,vCnt,l) v+(l>=vCnt{1})*(v<vCnt{2})*(mod(l,vCnt{3})==0)*vCnt{4};  % apply continuation
% ------------------------------------------------- PRE. 5) PREPARE FILTER
[dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
h = max( 0, rmin - sqrt( dx.^2 + dy.^2 ) );                                % conv. kernel
Hs = imfilter( ones( nely, nelx ), h, bcF );                               % matrix of weights (filter)
dHs = Hs;
% ------------------------ PRE. 6) ALLOCATE AND INITIALIZE OTHER PARAMETERS
[ x, dsK, dV ] = deal( zeros( nEl, 1 ) );                                  % initialize vectors
dV( act, 1 ) = 1/nEl/volfrac;                                              % derivative of volume (constant)
x( act ) = ( volfrac*( nEl - length(pasV) ) - length(pasS) )/length( act );% volume fraction on active set
x( pasS ) = 1;                                                             % set x = 1 on pasS set
[ xPhys, xOld, ch, loop, U ] = deal( x, 1, 1, 0, zeros( nDof, 1 ) );       % old x, x change, it. counter, U
% ================================================= START OPTIMIZATION LOOP
while ch > 1e-6 && loop < maxit
  loop = loop + 1;                                                         % update iter. counter
  % ----------- RL. 1) COMPUTE PHYSICAL DENSITY FIELD (AND ETA IF PROJECT.)
  xTilde = imfilter( reshape( x, nely, nelx ), h, bcF ) ./ Hs;
  xPhys( act ) = xTilde( act );                                            % reshape to column vector
  if ft > 1                              % compute optimal eta* with Newton
      f = ( mean( prj( xPhys, eta, beta ) ) - volfrac ) * ( ft == 3 );     % function (volume)
      while abs( f ) > 1e-6           % Newton process for finding opt. eta
          eta = eta - f / mean( deta( xPhys( : ), eta, beta ) );
          f = mean( prj( xPhys, eta, beta ) ) - volfrac;
      end
      dHs = Hs ./ reshape( dprj( xTilde, eta, beta ), nely, nelx );        % modification of the sensitivity
      xPhys = prj( xPhys, eta, beta );                                     % projected (physical) field
  end
  ch = norm( xPhys - xOld ) ./ sqrt( nEl );
  xOld = xPhys;
  % -------------------------- RL. 2) SETUP AND SOLVE EQUILIBRIUM EQUATIONS
  sK = ( Emin + xPhys.^penal * ( E0 - Emin ) );                            % stiffness interpolation
  dsK( act ) = -penal * ( E0 - Emin ) * xPhys( act ) .^ ( penal - 1 );     % derivative of stiffness interp.
  sK = reshape( Ke( : ) * sK', length( Ke ) * nEl, 1 );
  K = fsparse( Iar( :, 1 ), Iar( :, 2 ), sK, [ nDof, nDof ] );             % assemble stiffness matrix
  U( free ) = decomposition( K( free, free ), 'chol','lower' ) \ F( free );% solve equilibrium system
  % ------------------------------------------ RL. 3) COMPUTE SENSITIVITIES
  dc = dsK .* sum( ( U( cMat ) * Ke0 ) .* U( cMat ), 2 );                  % derivative of compliance
  dc = imfilter( reshape( dc, nely, nelx ) ./ dHs, h, bcF );               % filter objective sensitivity
  dV0 = imfilter( reshape( dV, nely, nelx ) ./ dHs, h, bcF );              % filter compliance sensitivity
  % ----------------- RL. 4) UPDATE DESIGN VARIABLES AND APPLY CONTINUATION
  xT = x( act );
  [ xU, xL ] = deal( xT + move, xT - move );                               % current upper and lower bound
  ocP = xT .* real( sqrt( -dc( act ) ./ dV0( act ) ) );                    % constant part in resizing rule
  l = [ 0, mean( ocP ) / volfrac ];                                        % initial estimate for LM
  while ( l( 2 ) - l( 1 ) ) / ( l( 2 ) + l( 1 ) ) > 1e-4                   % OC resizing rule
      lmid = 0.5 * ( l( 1 ) + l( 2 ) );
      x( act ) = max( max( min( min( ocP / lmid, xU ), 1 ), xL ), 0 );
      if mean( x ) > volfrac, l( 1 ) = lmid; else, l( 2 ) = lmid; end
  end
  [penal,beta] = deal(cnt(penal,penalCnt,loop), cnt(beta,betaCnt,loop));   % apply conitnuation on parameters
  % -------------------------- RL. 5) PRINT CURRENT RESULTS AND PLOT DESIGN
  fprintf( 'It.:%5i C:%7.4f V:%7.3f ch.:%0.2e penal:%7.2f beta:%7.1f eta:%7.2f \n', ...
      loop, F'*U, mean( xPhys ), ch, penal, beta, eta ); 
  colormap( gray ); imagesc( 1 - reshape( xPhys, nely, nelx ) );
  caxis([0 1]); axis equal off; drawnow;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by F. Ferrari, O. Sigmund                   %
% Dept. of Solid Mechanics-Technical University of Denmark,2800 Lyngby (DK)%
% Please send your comments to: feferr@mek.dtu.dk                          %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper Ferrari, F. Sigmund, O. - A new generation 99 %
% line Matlab code for compliance Topology Optimization and its extension  %
% to 3D, SMO, 2020                                                         %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - BOUNDARY CONDITIONS AND PASSIVE DOMAINS FOR THE TRUSS REINFORCEMENT EX.
% elNrs = reshape( 1 : nEl, nely, nelx );
% [lDofv,lDofh]=deal(2*nodeNrs(1,:),2*nodeNrs(:,end)-1);
% fixed = [1,2,nDof];
% a1=elNrs(1:nely/50,:);
% a2=elNrs( :,[1:nelx/50,end-nelx/50+1:end]);
% a3=elNrs(2*nely/5:end,2*nelx/5:end-nelx/5);
% [pasS,pasV]=deal(unique([a1(:);a2(:)]),a3(:));
%F = fsparse(lDofv',1,-2/(nelx+1),[nDof,1]) + ...
%    fsparse(lDofh,1,[0:1/(nely).^2:1/nely]',[nDof,1]);
