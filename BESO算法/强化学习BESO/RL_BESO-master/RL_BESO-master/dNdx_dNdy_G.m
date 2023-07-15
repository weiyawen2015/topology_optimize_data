function [dshapedx,dshapedy,G] = dNdx_dNdy_G(invjacob,dshapedxi,dshapedeta,eleNum)
%--------------------------------------------------------------------------
%  Computation of derivatives matrix G for stress matrix S for tangent
%  stiffness matrix K_T
%--------------------------------------------------------------------------

dshapedx = zeros(eleNum,4*4);
dshapedy = zeros(eleNum,4*4);
G = zeros(eleNum,4*4*8);
% G0 = zeros(4,8);
for i = 1:eleNum
    for j = 1:4
        temp = invjacob(eleNum,(j-1)*4+1:j*4);
        invjacob1 = reshape(temp,2,2);    
        
        dNdxi = dshapedxi(j,:);
        dNdeta = dshapedeta(j,:);
        
        temp = invjacob1*[dNdxi;dNdeta];
        dNdx = temp(1,:); 
        dNdy = temp(2,:);
        
        G1 = zeros(4,8);
        G1(1,1:2:end-1) = dNdx;
        G1(2,1:2:end-1) = dNdy;
        G1(3,2:2:end) = dNdx;
        G1(4,2:2:end) = dNdy;
        
        dshapedx(i,(j-1)*4+1:j*4) = dNdx;
        dshapedy(i,(j-1)*4+1:j*4) = dNdy;
        G(i,(j-1)*4*8+1:j*4*8) = G1(:)';
    end       
end
end