%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%


clear 
clc
nelx=60;   % x�᷽��ĵ�Ԫ��
nely=20;   % y�᷽��Ԫ��
volfrac=0.5;  %�����
penal=3;      %���ϲ�ֵ�ĳͷ�����
rmin=2.4;       %���ȹ��˵İ뾶
ft = 2;




%  function top88(nelx,nely,volfrac,penal,rmin,ft)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);%��Ԫ�ն���
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);%(nely+1)*(nelx+1)���ڵ�
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);%ȡ��
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);%ȡ��
%% PREPARE FILTER   �˲�
  % Դ����
        iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
        jH = ones(size(iH));
        sH = zeros(size(iH));
        k = 0;
        for i1 = 1:nelx
          for j1 = 1:nely
            e1 = (i1-1)*nely+j1;
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
              for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
              end
            end
          end
        end
        H = sparse(iH,jH,sH);
        Hs = sum(H,2);% �ܸ�  ϡ�����
  % ���
%     [dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1, ...
%     -ceil(rmin)+1:ceil(rmin)-1);
%     h = max(0,rmin-sqrt(dx.^2+dy.^2));
%     Hs = conv2(ones(nely,nelx),h,'same');
  %PDE
%         Rmin=rmin/2/sqrt(3);
%         KEF=Rmin^2*[4 -1 -2 -1; -1 4 -1 -2;...
%         -2 -1 4 -1; -1 -2 -1 4]/6 + ...
%         [4 2 1 2; 2 4 2 1;...
%         1 2 4 2; 2 1 2 4]/36;
%         edofVecF=reshape(nodenrs(1:end-1,1:end-1),nelx*nely,1);
%         edofMatF=repmat(edofVecF,1,4)...
%         +repmat([0 nely+[1:2] 1],nelx*nely,1);
%         iKF=reshape(kron(edofMatF,ones(4,1))',16*nelx*nely,1);
%         jKF=reshape(kron(edofMatF,ones(1,4))',16*nelx*nely,1);
%         sKF=reshape(KEF(:)*ones(1,nelx*nely),16*nelx*nely,1);
%         KF=sparse(iKF,jKF,sKF);
%         LF=chol(KF,'lower');
%         iTF=reshape(edofMatF,4*nelx*nely,1);
%         jTF=reshape(repmat([1:nelx*nely],4,1)',4*nelx*nely,1);
%         sTF=repmat(1/4,4*nelx*nely,1);
%         TF=sparse(iTF,jTF,sTF);
%% INITIALIZE ITERATION iterationѭ�� ����
x = repmat(volfrac,nely,nelx);%������ʼֵ  ����
xPhys = x;
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);%������Ԫ����
  %% objective function and sensitivity analysisĿ�꺯�� �������ȷ���
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)); % Ŀ�꺯��
  dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;%�����ȷ���
  dv = ones(nely,nelx);%V����Ʊ�����������һֱ��1
  %% FILTERING/MODIFICATION OF SENSITIVITIES    ����������sensitivity
  %Դ����
          if ft == 1
            dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
          elseif ft == 2
            dc(:) = H*(dc(:)./Hs);
            dv(:) = H*(dv(:)./Hs);
          end
   %���
%          if ft == 1
%             dc = conv2(dc.*xPhys,h,'same')./Hs./max(1e-3,xPhys);
%             elseif ft == 2
%             dc = conv2(dc./Hs,h,'same');
%             dv = conv2(dv./Hs,h,'same');
%          end
    %PDE
%             if ft == 1
%             dc(:) = (TF'*(LF'\(LF\(TF*(dc(:).*xPhys(:))))))...
%             ./max(1e-3,xPhys(:));
%             elseif ft == 2
%             dc(:) = TF'*(LF'\(LF\(TF*dc(:))));
%             dv(:) = TF'*(LF'\(LF\(TF*dv(:))));
%             end
  %% OC׼�� OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES criterion��׼variables���� 
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    
    
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
       xPhys(:) = (H*xnew(:))./Hs;                  %Դ����
%         xPhys = conv2(xnew,h,'same')./Hs;            %���
%         xPhys(:) = (TF'*(LF'\(LF\(TF*xnew(:)))));    %PDE
    end
    
    % ���ַ�
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end

% 
% clear 
% clc
% nelx=60;   % x�᷽��ĵ�Ԫ��
% nely=20;   % y�᷽��Ԫ��
% volfrac=0.5;  %�����
% penal=3;      %���ϲ�ֵ�ĳͷ�����
% rmin=1.5;       %���ȹ��˵İ뾶


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

