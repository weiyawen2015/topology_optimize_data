%%%% A SOFT-KILL BESO CODE BY X. HUANG and Y.M. XIE %%%%
% function besoall(nelx,nely,volfrac,er,rmin)
clear;
clc;
close all
format short;

nelx = 80;
nely = 40
volfrac = 0.3;
er = 0.01;
rmin = 2.5;



% INITIALIZE
x(1:nely,1:nelx) = 1.0;
vol=1.0; 
i = 0; 
change = 1.0;
penal = 3.0; 

% START iTH ITERATION
while change > 0.0001
% for i = 0:5
    i
    
    
i = i + 1; 
vol = max(vol*(1-er),volfrac);
if i >1;
    olddc = dc; 
end
% FE-ANALYSIS
[U]=FE(nelx,nely,x,penal);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
[KE] = lk;
c(i) = 0.0;
for ely = 1:nely
    for elx = 1:nelx
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)*elx+ely;
        Ue = U([2*n1-1;2*n1;2*n2-1;2*n2;2*n2+1;2*n2+2;2*n1+1;2*n1+2],1);
        c(i) = c(i) + 0.5*x(ely,elx)^penal*Ue'*KE*Ue;
        dc(ely,elx) = 0.5*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
    end
end
% FILTERING OF SENSITIVITIES
[dc] = check(nelx,nely,rmin,x,dc);
% STABLIZATION OF EVOLUTIONARY PROCESS
if i > 1;
     dc = (dc+olddc)/2; 
end

% BESO DESIGN UPDATE
[x] = ADDDEL(nelx,nely,vol,dc,x);
% PRINT RESULTS
if i>10;
    change=abs(sum(c(i-9:i-5))-sum(c(i-4:i)))/sum(c(i-4:i));
end
disp([ ' It.: ' sprintf('%4i',i) ' Obj.: ' sprintf('%10.4f',c(i)) ...
    ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
    ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES
colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
end



