function p = DIRT_PModel_Compose(pmodel, p2, p1, color)
% DIRT_PModel_Compose(pmodel, p2, p1, color)
% 
% p = p2 \circ p1 
% p(v) = p2(p1(v))
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~color
    switch pmodel
        case 'GainAndBias'
            p = [ 
                p2(1)*p1(1)
                p2(1)*p1(2)+p2(2)
                ];
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
else
    switch pmodel
        case 'SingleGainAndBias'
            p = [ 
                p2(1)*p1(1)
                p2(1)*p1(2)+p2(2)
                ];
        case 'MultipleGainsAndBiases'
            p = [ 
                p2(1)*p1(1)
                p2(2)*p1(2)
                p2(3)*p1(3)
                p2(1)*p1(4)+p2(4)
                p2(2)*p1(5)+p2(5)
                p2(3)*p1(6)+p2(6)
                ];
        case 'Affine'
            A1 = reshape(p1(1:9), 3, 3);
            A2 = reshape(p2(1:9), 3, 3);
            A = A2 * A1;
            b = A2 * p1(10:12) + p2(10:12);
            p = [ A(:) ; b ];
        otherwise
            error(['unknown photometric model ' pmodel]);
    end;
end;


