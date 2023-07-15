%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%
function [x]=ADDDEL(nelx,nely,vol,dc,x)
l1 = min(min(dc)); l2 = max(max(dc));

   figure(1)
   subplot(1,2,1)
   colormap(gray); imagesc(-dc); axis equal; axis tight; axis off;
%     colorbar;

while ((l2-l1)/l2 > 1.0e-2)
    th = (l1+l2)/2.0;
    x = max(0.001,sign(dc-th));
   
    subplot(1,2,2)
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off;colorbar;
    
    if sum(sum(x))-vol*(nelx*nely) > 0;
        l1 = th;
    else
        l2 = th;
    end
    title( ['vol=',num2str(roundn(vol,-3)),'VOL=',num2str(roundn(sum(sum(x))/(nelx*nely),-3)),';l1=',num2str(roundn(l1,-3)),';l2=',num2str(roundn(l2,-3))])
end