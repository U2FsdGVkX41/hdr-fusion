function p = DIRT_PModel_Identity(pmodel, color)
% p = DIRT_PModel_Identity(pmodel, color)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~color
    switch pmodel
        case 'GainAndBias'
            p = [ 1 0 ]';
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
else
    switch pmodel
        case 'SingleGainAndBias'
            p = [ 1 0 ]';
        case 'MultipleGainsAndBiases'
            p = [ 1 1 1 0 0 0 ]';
        case 'Affine'
            p = [ 1 0 0 0 1 0 0 0 1 0 0 0 ]';
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
end;
