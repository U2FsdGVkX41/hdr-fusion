function DIRT_PModel_DispInit(pmodel, color)
% DIRT_PModel_DispInit(pmodel, color)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~color
    switch pmodel
        case 'GainAndBias'
            fprintf(' | gain         | bias');
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
else
    switch pmodel
        case 'SingleGainAndBias'
            fprintf(' | gain         | bias');
        case 'MultipleGainsAndBiases'
            fprintf(' | mean gain    | mean bias');
        case 'Affine'
            fprintf(' | mean gain    | mean bias');
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
end;


