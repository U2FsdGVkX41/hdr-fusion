function p = DIRT_PModel_Update(pmodel, color, p, dp)
% p = DIRT_PModel_Update(pmodel, color, p, dp)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~color
    switch pmodel
        case 'GainAndBias'
            invda = 1 / ( 1 + dp(1) );
            p(1) = invda * p(1);
            p(2) = invda * (p(2) - dp(2));
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
else
    switch pmodel
        case 'SingleGainAndBias'
            invda = 1 / ( 1 + dp(1) );
            p(1) = invda * p(1);
            p(2) = invda * (p(2) - dp(2));
        case 'MultipleGainsAndBiases'
            invda = 1 ./ ( 1 + dp(1:3) );
            p(1:3) = invda .* p(1:3);
            p(4:6) = invda .* ( p(4:6) - dp(4:6) );
        case 'Affine'
            invdA = inv( eye(3) + reshape(dp(1:9), 3, 3) );
            p(1:9) = vect( invdA * reshape(p(1:9), 3, 3) );
            p(10:12) = invdA * ( p(10:12) - dp(10:12) );
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
end;


