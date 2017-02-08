function np = DIRT_PModel_N(pmodel, color)
% np = DIRT_PModel_N(pmodel, color)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~color
    switch pmodel
        case 'GainAndBias'
            np = 2;
        otherwise
            error(['unknown photometric model ' pmodel]);
    end;
else
    switch pmodel
        case 'SingleGainAndBias'
            np = 2;
        case 'MultipleGainsAndBiases'
            np = 6;
        case 'Affine'
            np = 12;
        otherwise
            error(['unknown photometric model ' pmodel]);
    end;
end;


