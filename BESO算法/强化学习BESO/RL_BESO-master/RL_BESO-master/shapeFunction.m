function [wt,shape,dshapedxi,dshapedeta] = shapeFunction()
%------------------------------------------------------------------------
%  Purpose:
%     compute isoparametric four-node Quadilateral shape functions
%     and their derivatves at the selected (integration) point
%     in terms of the natural coordinate 

%
%  Variable Description:
%     shape - shape functions for four-node element
%     dshapedxi - derivatives of the shape functions w.r.t. xi
%     dshapedeta - derivatives of the shape functions w.r.t. eta
%
%  Notes:
%     1st node at (-1,-1), 2nd node at (1,-1)
%     3rd node at (1,1), 4th node at (-1,1)
%------------------------------------------------------------------------
gausspoint = [-0.577350269189626 -0.577350269189626;
    0.577350269189626  -0.577350269189626;
    0.577350269189626 0.577350269189626;
    -0.577350269189626  0.577350269189626] ;
wt = [1 ;1;1; 1];


xi = gausspoint(:,1);
eta = gausspoint(:,2);

% shape functions, [N1,N2,N3,N4;] for every gauss point
shape = zeros(4,4);
shape(:,1) = 0.25.*(1-xi).*(1-eta);
shape(:,2) = 0.25.*(1+xi).*(1-eta);
shape(:,3) = 0.25.*(1+xi).*(1+eta);
shape(:,4) = 0.25.*(1-xi).*(1+eta);

% derivatives, [dN1dxi,dN2dxi,dN3dxi,dN4dxi;] for every gauss point
dshapedxi = zeros(4,4);
dshapedxi(:,1)=-0.25*(1-eta);
dshapedxi(:,2)=0.25*(1-eta);
dshapedxi(:,3)=0.25*(1+eta);
dshapedxi(:,4)=-0.25*(1+eta);

%[dN1deta,dN2deta,dN3deta,dN4deta;] for every gauss point
dshapedeta = zeros(4,4);
dshapedeta(:,1)=-0.25*(1-xi);
dshapedeta(:,2)=-0.25*(1+xi);
dshapedeta(:,3)=0.25*(1+xi);
dshapedeta(:,4)=0.25*(1-xi);
end