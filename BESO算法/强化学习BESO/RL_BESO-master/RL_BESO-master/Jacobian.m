function [detjacob,invjacob] = Jacobian(nodes,eles,dshapedxi,dshapedeta,eleNum)
%--------------------------------------------------------------------------
% for gauss integration
%--------------------------------------------------------------------------

detjacob = zeros(eleNum,4);		% area
invjacob = zeros(eleNum,4*4);	
for i = 1:eleNum
    node = eles(i,:);
    coord = nodes(node,:);
    
    % for 4 gauss points
    for j = 1:4
        dNdxi = dshapedxi(j,:);
        dNdeta = dshapedeta(j,:);
    
        jacobian = [sum(dNdxi'.*coord(:,1)),sum(dNdxi'.*coord(:,2));...
                    sum(dNdeta'.*coord(:,1)),sum(dNdeta'.*coord(:,2))];
            
        detjacobian = det(jacobian) ;  % Determinant of Jacobian matrix
        invjacobian = inv(jacobian) ;  % Inverse of Jacobian matrix
        detjacob(i,j) = detjacobian;
        invjacob(i,(j-1)*4+1:j*4) = invjacobian(:)';
%         jacob(eleNum,(j-1)*5+1:j*5) = [detjacobian,invjacobian(:)'];
    end
end
end