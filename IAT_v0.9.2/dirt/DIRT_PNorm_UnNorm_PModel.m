function p = DIRT_PNorm_UnNorm_PModel(params, pmodel, color, p)
% p = DIRT_PNorm_UnNorm_PModel(params, pmodel, color, p)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~color
    switch pmodel
        case 'GainAndBias'
            p(2) = p(2) * params;
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
else
    switch pmodel
        case 'SingleGainAndBias'
            p(2) = p(2) * params;
        case 'MultipleGainsAndBiases'
            p(4:6) = p(4:6) * params;
        case 'Affine'
            p(10:12) = p(10:12) * params;
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
end;


