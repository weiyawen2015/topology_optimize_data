%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc)  
L1 = 0; L2 = 100000; move = 0.2;
while (L2-L1 > 1e-4) %l 是小写的L 不是12 11  
    
  Lmid = 0.5*(L2+L1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./Lmid)))));
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;
    L1 = Lmid;
  else
    L2 = Lmid;
  end
end