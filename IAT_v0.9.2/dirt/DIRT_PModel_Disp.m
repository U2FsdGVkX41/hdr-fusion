function DIRT_PModel_Disp(pmodel, color, p)
% DIRT_PModel_Disp(pmodel, color, p)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~color
    switch pmodel
        case 'GainAndBias'
            fprintf(' | %012.7f | %012.7f', p(1), p(2));
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
else
    switch pmodel
        case 'SingleGainAndBias'
            fprintf(' | %012.7f | %012.7f', p(1), p(2));
        case 'MultipleGainsAndBiases'
            fprintf(' | %012.7f | %012.7f', mean(abs(p(1:3))), mean(abs(p(4:6))));            
        case 'Affine'
            fprintf(' | %012.7f | %012.7f', mean(abs(p(1:9))), mean(abs(p(10:12))));                        
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
end;


