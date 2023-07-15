clear all
close all
clc
nelx=60;   % x�᷽��ĵ�Ԫ��
nely=20;   % y�᷽��Ԫ��
volfrac=0.5;  %�����
% penal=3;      %���ϲ�ֵ�ĳͷ�����
rmin=2.4;       %���ȹ��˵İ뾶
alfa=2.5;
xfar = 1;




xshi =0:0.01:1;
% flag==1
  tmpfx=(xshi-xfar);
 fx=(1+tanh(alfa*tmpfx));
 fx2 = xshi.^3;
  figure(1)
  plot(xshi,fx,'b*')
  hold on
   plot(xshi,fx2,'ro') 
  hold off






% INITIALIZE
x(1:nely,1:nelx) = volfrac; % x����Ʊ��� ����ʼֵ
loop = 0; 
change = 1.;
% START ITERATION
 while change > 0.01  

  loop = loop + 1;
  xold = x; % �洢ǰһ�α���
  
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,alfa,xfar);      %��������Ԫ����      
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  [KE] = lk;  %��Ԫ�նȾ���
  
  c = 0.; %�������Ŀ�꺯���ı���.����Ŀ�꺯���Ǹն����Ҳ������
                                     %����С  
   % ��Ԫ����ɨ��                                  
  for ely = 1:nely  
    for elx = 1:nelx
      n1 = (nely+1)*(elx-1)+ely; 
      n2 = (nely+1)* elx   +ely;
      Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);% 8x1�ľ���Ԫ��λ������
%       c = c + x(ely,elx)^penal*Ue'*KE*Ue;       %����Ŀ�꺯����ֵ��������Ƕ��٣�
% 
        c = c + (1+tanh(alfa*( x(ely,elx)-xfar)))*Ue'*KE*Ue;       %����Ŀ�꺯����ֵ��������Ƕ���
        dc(ely,elx) = -alfa*(1-tanh(alfa*( x(ely,elx)-xfar))^2)*Ue'*KE*Ue; % �����ȷ����Ľ�� Ϊ��Ʊ���������׼��
   
        
%         c = c + (1/(1+exp(-alfa*(x(ely,elx)-0.5))))*Ue'*KE*Ue;       %����Ŀ�꺯����ֵ��������Ƕ���
%         dc(ely,elx) = -alfa*exp(-alfa*(x(ely,elx)-0.5))/(1+exp(-alfa*(x(ely,elx)-0.5)))^2*Ue'*KE*Ue; % �����ȷ����Ľ�� Ϊ��Ʊ���������׼��     
    
    end
  end
  
  
% FILTERING OF SENSITIVITIES
   [dc]   = check(nelx,nely,rmin,x,dc);          %�����ȹ��ˣ�Ϊ�˱߽��˳һ��
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [x]    = OC(nelx,nely,x,volfrac,dc);       %�Ż�׼�򷨸�����Ʊ��� ������xnew
  save x x
% PRINT RESULTS
  change = max(max(abs(x-xold))); %�µ���Ʊ�����ɵ���Ʊ����Ĳ�ľ���ֵ�����ֵ
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])          %��Ļ����ʾ������Ϣ
% PLOT DENSITIES  
  figure(2)
  colormap(gray); % x�Ǿ��� ��������ͼ�� ��"-"��Ϊ�˺ڰ׵ߵ�һ��
  imagesc(-x);  %colormap(gray);imagesc(-x);���ʹ�� ��ϸ���Ϳ�https://blog.csdn.net/zy122121cs/article/details/49761307
  axis equal; 
  axis tight; % ��ʾ����һЩ
  axis off;
  pause(1e-7); 
%    title([' ��loop��ѭ��' num2str(loop)]) 
  %�Ż������ͼ����ʾ��������Ϊ����ͼ����ʾ�����ܲ��ã�̫���ˡ��ȽϷ����ͼ����ʾӦ���ǣ�
% ÿһ�ε���ͬʱ��ʾ�Ż������Ŀ�꺯�����ߣ�Ȼ���Զ�����ÿһ�εĽ����

end
