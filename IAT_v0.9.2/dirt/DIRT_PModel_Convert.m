function pTo = DIRT_PModel_Convert(pmodelFrom, colorFrom, pmodelTo, colorTo, pFrom)
% DIRT_PModel_Convert(pmodelFrom, colorFrom, pmodelTo, colorTo, pFrom)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~colorFrom
    switch pmodelFrom
%        case 'GainAndBias'
%            p = [ 1 0 ]';
        otherwise
            error(['unknown photometric model' pmodelFrom]);
    end;
else
    switch pmodelFrom
%        case 'SingleGainAndBias'
%            p = [ 1 0 ]';
        case 'MultipleGainsAndBiases'
            if ~colorTo
                switch pmodelTo
                    otherwise
                        error(['unknown photometric model' pmodelTo]);
                end;
            else
                switch pmodelTo
                    case 'Affine'
                        pTo = [ pFrom(1) 0 0 0 pFrom(2) 0 0 0 pFrom(3) pFrom(4:6)' ]';
                    otherwise
                        error(['unknown photometric model' pmodelTo]);
                end;
            end;
%        case 'Affine'
%            p = [ 1 0 0 0 1 0 0 0 1 0 0 0 ]';
        otherwise
            error(['unknown photometric model' pmodelFrom]);
    end;
end;
