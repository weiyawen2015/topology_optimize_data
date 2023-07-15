function [U]=FE(nelx,nely,x,alfa,xfar)
[KE] = lk; %��Ԫ�նȾ���
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
    %����նȾ����ϡ����� ע��ڵ����͵�Ԫ���Ĺ�ϵ
F = sparse(2*(nely+1)*(nelx+1),1); %�غ����� ����ô��ڵ�(nely+1)*(nelx+1)
U = zeros(2*(nely+1)*(nelx+1),1); %λ������ ����ô��ڵ�(nely+1)*(nelx+1)
% �ܸ���װ
for elx = 1:nelx   % ��elx�еĵ�Ԫ
  for ely = 1:nely  %��ely�еĵ�Ԫ
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];%�����Y���Ƿ���ģ����ǲ�Ӱ�����Ľ�������������¥TYNGOD��λ���ֵĽ��ͣ���лTYNGOD��
%     K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;  %����Ԫ�նȾ�����װ���ܵĸնȾ���
       K(edof,edof) = K(edof,edof) +  (1+tanh(alfa*( x(ely,elx)-xfar)))*KE;  %����Ԫ�նȾ�����װ���ܵĸնȾ���
%             K(edof,edof) = K(edof,edof) + (1/(1+exp(-alfa*(x(ely,elx)-0.5))))*KE;  %����Ԫ�նȾ�����װ���ܵĸնȾ���
  
  
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F(2,1) = -1; % ����Ԫ�ض���0

fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]); %unionȡ������̶����
alldofs     = [1:2*(nely+1)*(nelx+1)];                       %���н�� 
freedofs    = setdiff(alldofs,fixeddofs);  % setdiff����A���У�B��û�е�ֵ���ɽڵ�
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:); %���   
U(fixeddofs,:)= 0;












