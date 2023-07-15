%%%% CODE of Combining BESO with RL BY H.B.SUN and L.MA in 2019 %%%%
function [ITER,C,IOU,C_difference,nsubopt]=RESO_cantilever(nelx,nely,volfrac,er,rmin,sym,emax,w);
%% INITIALIZE
x(1:nely,1:nelx) = 1.; vol=1.; i = 0; change = 1.; penal = 3.;Emin=10^(-9);Emax=1.0;
pic=zeros(100,1);R=zeros(100,nelx*nely);V=zeros(100,nelx*nely);gamma=0.9;winx=1;winy=1;
Nsubopt=20;unfixvol=1-volfrac/2;
C=zeros(Nsubopt+1,1);ITER=zeros(Nsubopt+1,1);C_difference=zeros(Nsubopt,1);IOU=zeros(Nsubopt,1);
%% PREPARE FILTER
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
Hs = sum(H,2);
%% START iTH ITERATION
while (change > 0.001)&&(i<100)
    i = i + 1; 
    vol = max(vol*(1-er),volfrac);
    if i >1; 
        olddc = dc; 
    end
% FE-ANALYSIS
    [U]=FE(nelx,nely,x,penal);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [KE] = lk;
    c(i) = 0.;
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1)+ely;
            n2 = (nely+1)* elx +ely;
            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2;2*n1+1;2*n1+2],1);
%             c(i) = c(i) + 0.5*x(ely,elx)^penal*Ue'*KE*Ue;
%             dc(ely,elx) = 0.5*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
            c(i) = c(i) + 0.5*(Emin+x(ely,elx)^penal*(Emax-Emin))*Ue'*KE*Ue;
            dc(ely,elx) = 0.5*x(ely,elx)^(penal-1)*(Emax-Emin)*Ue'*KE*Ue;            
        end
    end
% FILTERING OF SENSITIVITIES
    dc(:) = H*(dc(:)./Hs);
%     dv(:) = H*(dv(:)./Hs);
%     [dc] = check(nelx,nely,rmin,x,dc);
% STABLIZATION OF EVOLUTIONARY PROCESS
    if i > 1; 
        dc = (dc+olddc)/2.; 
    end
% BESO DESIGN UPDATE
    [x] = ADDDEL(nelx,nely,vol,dc,x);
    if (vol-volfrac)>0.01
        qi=i+5;
    end    
    if qi>=i
        R(i,:)=reshape(x,1,nelx*nely);
    end
% PRINT RESULTS
    if i>10
        change=abs(sum(c(i-9:i-5))-sum(c(i-4:i)))/sum(c(i-4:i));
    end
    volnext(i)=mean(x(:));
    disp([' It.: ' sprintf('%4i',i)  ' Obj.: ' sprintf('%10.4f',c(i)) ' Vol.: ' sprintf('%6.3f',volnext(i)) ' ch.: ' sprintf('%6.3f',change )])
% % PLOT DENSITIES
    figure(1);
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off; 
end
% % plot vol    
% figure(60);
% yyaxis right;
% plot(1:i,volnext(1:i),'*'); 
% xlabel('Iteration')
% ylabel('Volume Fraction')
% plot dc    
figure(60);hold on;
plot(1:i,c(1:i),'color','b','LineWidth',2);  
xlabel('Iteration')
ylabel('Compliance')
c0=c(i);
C(1)=c0;ITER(1)=i;
xsubopt=reshape(x,1,nelx*nely);nsubopt=[];N=0;
%% RL section
for z=1:Nsubopt
    if c(i)/c0-1<0.2
        N=N+1;
        for j=qi-1:-1:1
            R(j,:)=R(j,:)+gamma*R(j+1,:);
        end
        V(1:qi,:)=V(1:qi,:)+1/N*(R(1:qi,:)*(2-c(i)/c0)-V(1:qi,:));
    end
%     % plot V curve
%     for i=1:5:16
%         [V0,index]=sort(V(i,:));
%         figure(20);hold on;
%         if i==1
%             plot(1:nely*nelx,V0);
%         elseif i==6
%             plot(1:nely*nelx,V0,':');
%         elseif i==11
%             plot(1:nely*nelx,V0,'--');
%         else
%             plot(1:nely*nelx,V0,'-.');
%         end
%     end
%     legend('1st iteration','6th iteration','11th iteration','16th iteration');
%     xlabel('number of elements');
%     ylabel('Value Function');
% INITIALIZE
x = ones(nely,nelx); vol=1.; i = 0; change = 1.; 
R=zeros(100,nelx*nely);
% START iTH ITERATION
while (change > 0.001)&&(i<100)||(vol>volfrac)
    i = i + 1; 
    vol = max(vol*(1-er),volfrac);
    if i >1; 
        olddc = dc; 
    end
% FE-ANALYSIS
    [U]=FE(nelx,nely,x,penal);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [KE] = lk;
    c(i) = 0.;
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1)+ely;
            n2 = (nely+1)* elx +ely;
            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2;2*n1+1;2*n1+2],1);
            c(i) = c(i) + 0.5*x(ely,elx)^penal*Ue'*KE*Ue;
            dc(ely,elx) = 0.5*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
        end
    end
% FILTERING OF SENSITIVITIES
    dc(:) = H*(dc(:)./Hs);
% STABLIZATION OF EVOLUTIONARY PROCESS
    if i > 1; 
        dc = (dc+olddc)/2.; 
    end
% RL or BESO DESIGN UPDATE
    e=emax+(0-emax)*min(1,1/(qi-5)*i)^3;
    if (vol-volfrac)>0.01;
%     if (vol-volfrac)>0.01||qi>=i;
        [x] = RL_ADDDEL(i,qi,nelx,nely,unfixvol,volnext(i),V(i,:),dc,e,winx,winy,sym,w);
    else
        [x] = ADDDEL(nelx,nely,vol,dc,x);
    end
    if qi>=i
        R(i,:)=reshape(x,1,nelx*nely);
    end
% PRINT RESULTS
    if i>10;
        change=abs(sum(c(i-9:i-5))-sum(c(i-4:i)))/sum(c(i-4:i));
    end
    disp([' It.: ' sprintf('%4i',i)  ' Obj.: ' sprintf('%10.4f',c(i)) ' Vol.: ' sprintf('%6.3f',mean(x(:))) ' ch.: ' sprintf('%6.3f',change )])

%     % plot dc curve
%     if i==1||i==6||i==11||i==16
%         dc1=reshape(dc,1,nelx*nely);
%         [dc0,index]=sort(dc1);
%         if i==1
%             figure(30);hold on;
%             plot(1:nely*nelx,dc0(1:nelx*nely));
%             figure(40);hold on;
%             dc0=min(V(i,:))+(dc0-dc0(1))/(dc0(0.75*nelx*nely)-dc0(1))*(max(V(i,:))-min(V(i,:)));
%             plot(1:0.75*nely*nelx,dc0(1:0.75*nelx*nely));
%         elseif i==6
%             figure(30);hold on;
%             plot(1:nely*nelx,dc0(1:nelx*nely),':');
%             figure(40);hold on;
%             dc0=min(V(i,:))+(dc0-dc0(1))/(dc0(0.75*nelx*nely)-dc0(1))*(max(V(i,:))-min(V(i,:)));
%             plot(1:0.75*nely*nelx,dc0(1:0.75*nelx*nely),':');            
%         elseif i==11
%             figure(30);hold on;
%             plot(1:nely*nelx,dc0(1:nelx*nely),'--');
%             figure(40);hold on;
%             dc0=min(V(i,:))+(dc0-dc0(1))/(dc0(0.75*nelx*nely)-dc0(1))*(max(V(i,:))-min(V(i,:)));
%             plot(1:0.75*nely*nelx,dc0(1:0.75*nelx*nely),'--');    
%         else
%             figure(30);hold on;
%             plot(1:nely*nelx,dc0(1:nelx*nely),'-.');
%             figure(40);hold on;
%             dc0=min(V(i,:))+(dc0-dc0(1))/(dc0(0.75*nelx*nely)-dc0(1))*(max(V(i,:))-min(V(i,:)));
%             plot(1:0.75*nely*nelx,dc0(1:0.75*nelx*nely),'-.');    
%         end
%     end
end
% legend('1st iteration','6th iteration','11th iteration','16th iteration');
% xlabel('number of elements');
% ylabel('Sensitivity');
% figure(30);hold on;
% legend('1st iteration','6th iteration','11th iteration','16th iteration');
% xlabel('number of elements');
% ylabel('Sensitivity');
C(z+1)=c(i);ITER(z+1)=i;
C_difference(z)=C(z+1)/C(1)-1;
if C_difference(z)<0.2
    x=reshape(x,1,nelx*nely);
    for j=1:size(xsubopt,1)
        IOU(z)=max(IOU(z),length(find((x+xsubopt(j,:))>1.5))/length(find((x+xsubopt(j,:))>0.5)));
        if IOU(z)>0.9
            break;
        end
        if j==size(xsubopt,1)
            xsubopt=[xsubopt;x];
            nsubopt=[nsubopt;z+1];
            % plot c 
            if length(nsubopt)==1
                figure(60);hold on;
                plot(1:i,c(1:i),'k');
            elseif length(nsubopt)==3
                figure(60);hold on;
                plot(1:i,c(1:i),'r:');                
            elseif length(nsubopt)==5
                figure(60);hold on;
                plot(1:i,c(1:i),'g--');    
            elseif length(nsubopt)==8
                figure(60);hold on;
                plot(1:i,c(1:i),'m-.');    
                legend('BESO',['episode' sprintf('%4i',nsubopt(1))],['episode' sprintf('%4i',nsubopt(3))],['episode' sprintf('%4i',nsubopt(5))],['episode' sprintf('%4i',nsubopt(7))]);
                a=1;
            end
        end
    end
end
if length(nsubopt)==10
    break
end
end

for z=1:length(nsubopt)
    figure(1+z);
    x=reshape(xsubopt(1+z,:),nely,nelx);
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off;
end
% nindex=[];
% iter=1;
% for z=1:Nsubopt
%     if z==nsubopt(iter)
%         iter=iter+1;
%         nindex=[nindex;'accepted']
%     else
%         nindex=[nindex;'rejected']
%     end
% end

end

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=ADDDEL(nelx,nely,volfra,dc,x)
l1 = min(min(dc)); l2 = max(max(dc));
while ((l2-l1)/l2 > 1.0e-5)
    th = (l1+l2)/2.0;
    x = max(0.001,sign(dc-th));
    if mean(x(:))-volfra > 0;
        l1 = th;
    else
        l2 = th;
    end
end
end

%%%%%%%%%% RL UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=RL_ADDDEL(i,qi,nelx,nely,unfixvol,volnext,V,dc,e,winx,winy,sym,w)
% symmetric constraint
if sym==1
    V=reshape(V,nely,nelx);
    nely=nely/2;
    V=V(1:nely,:);
    V=reshape(V,1,nelx*nely);
    dc=dc(1:nely,:);
end 
% determine V by both V and dc
dc=reshape(dc,1,nelx*nely);
[dc0,index]=sort(dc);
Vnow=min(V)+(dc-dc0(1))/(dc0(unfixvol*nelx*nely)-dc0(1))*(max(V)-min(V));
V=w*V+(1-w)*Vnow;

% delete 1-e elements
[V0,index]=sort(V);
base=floor((1-e)*(1-volnext)*(nelx*nely));
x=ones(1,nely*nelx);
V_critical=V0(base); % RANDOM CHOICE
index_critical=find(V0==V_critical);
base2=index_critical(1);
index(index_critical)= index(index_critical(1)-1+randperm(length(index_critical)));
x(index(1:base))=0.001;
x=reshape(x,nely,nelx);
% delete e elements by window
% m=floor(e*(1-volnext)*(nelx*nely)/(2*win+1)^2);
m=floor(e*(1-volnext)*(nelx*nely)/(2*winx+1)/(2*winy+1));
if m>0
n=unfixvol*nelx*nely-base;
p=randperm(n);
for i=1:m
    ind=index(base+p(i));
    indx=ceil(ind/nely);
    indy=mod(ind,nely);
    if indy==0
        indy=nely;
    end
    x(indy,max(1,indx-winx):min(indx+winx,nelx))=0.001;
    if winy>0
        for j=1:winy
            x(max(1,indy-j),max(1,indx-winx):min(indx+winx,nelx))=0.001;
            x(min(nely,indy+j),max(1,indx-winx):min(indx+winx,nelx))=0.001;
        end
    end
end
end
% delete rest elements
n=floor((mean(x(:))-volnext)*(nelx*nely));
if n>0
    ind2=find(x(index(1:unfixvol*nelx*nely))==1);
    p=randperm(length(ind2));
    x(index(ind2(p(1:n))))=0.001;
end
if sym==1
    x=[x;flipud(x)];
end 
end

%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk;
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); 
U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
    for ely = 1:nely
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx +ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
    end
end
% DEFINE LOADS AND SUPPORTS (Cantilever)
fixeddofs=[1:2*(nely+1)]; %left
F = sparse(2*(nelx+1)*(nely+1)-nely,1,-1,2*(nely+1)*(nelx+1),1); %right middle
% F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1); %middle top
% F = sparse(2*(nely+1),1,-1,2*(nely+1)*(nelx+1),1); %middle bottom
% fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]); %left x and right bottom y
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs,:)= 0;
end

%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
%% MATERIAL PROPERTIES
E = 1;
nu = 0.3;
k=[1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8 nu/6 1/8-3*nu/8];

KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8); ...
k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3); ...
k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2); ...
k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5); ...
k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4); ...
k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7); ...
k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6); ...
k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end
