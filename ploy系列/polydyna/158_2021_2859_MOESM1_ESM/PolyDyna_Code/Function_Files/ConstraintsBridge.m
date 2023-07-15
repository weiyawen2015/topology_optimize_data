%-------------------------- ConstraintsBridge ----------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    %
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
VolFrac = 0.3;
ElemCtrd = Centroids(fem);
Hb = max(Node(:,2)); % Height of bridge
fem.SElem = find(ElemCtrd(:,2)>=Hb-L/40); % Elements in passive region
ElemInd{1} = setdiff((1:fem.NElem)',fem.SElem);
%-------------------------------------------------------------------------%