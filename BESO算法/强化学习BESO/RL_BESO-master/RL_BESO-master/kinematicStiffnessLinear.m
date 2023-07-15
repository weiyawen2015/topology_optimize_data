function B0 = kinematicStiffnessLinear(dshapedx,dshapedy,eleNum)
%--------------------------------------------------------------------------
% linear kinematic Stiffness matrix of strain
%--------------------------------------------------------------------------

B0 = zeros(eleNum,3*8*4);
for i = 1:eleNum
    for j = 1:4
        dNdx = dshapedx(i,(j-1)*4+1:j*4);
        dNdy = dshapedy(i,(j-1)*4+1:j*4);
        
        B0i = zeros(3,8);
        B0i(1,1:2:end-1) = dNdx;
        B0i(2,2:2:end) = dNdy;
        B0i(3,2:2:end) = dNdx;
        B0i(3,1:2:end-1) = dNdy;
        
        B0(i,(j-1)*3*8+1:j*3*8) = B0i(:)';
    end
end
end