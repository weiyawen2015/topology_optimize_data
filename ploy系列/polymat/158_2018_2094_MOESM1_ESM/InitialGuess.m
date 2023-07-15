%-------------------------------- PolyMat --------------------------------%
% Ref: ED Sanders, A Pereira, MA Aguilo, GH Paulino, "PolyMat: an         %
%      efficient Matlab code for multi-material topology optimization,"   %
%      Struct Multidisc Optim, DOI 10.1007/s00158-018-2094-0, 2018        %
%-------------------------------------------------------------------------%
function [zIni] = InitialGuess(VolFrac,ElemInd,MatInd,SElemInd,SMatInd,zIni)
NConstr = size(VolFrac,1);
for i = 1:NConstr
    MatIndices = cell2mat(MatInd(i));
    ElemIndices = cell2mat(ElemInd(i));
    for m = 1:size(MatIndices,2)
        zIni(ElemIndices,MatIndices(m)) = VolFrac(i)./size(MatIndices,2);
    end
end
NFixed = size(SElemInd,1);
for i = 1:NFixed
    MatIndices = cell2mat(SMatInd(i));
    ElemIndices = cell2mat(SElemInd(i));
    for m = 1:size(MatIndices,2)
        zIni(ElemIndices,MatIndices(m)) = 1;
    end
end
%-------------------------------------------------------------------------%