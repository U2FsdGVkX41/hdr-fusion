function S = DIRT_PModel_Apply(pmodel, color, p, S)
% DIRT_PModel_Apply(pmodel, color, p, S)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~color
    switch pmodel
        case 'GainAndBias'
            S = p(1) * S + p(2);
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
else
    switch pmodel
        case 'SingleGainAndBias'
            S = p(1) * S + p(2);
        case 'MultipleGainsAndBiases'
            % XXX do not work for 3D arrays
            if length(size(S)) == 2
                S(:, 1) = p(1) * S(:, 1) + p(4);
                S(:, 2) = p(2) * S(:, 2) + p(5);
                S(:, 3) = p(3) * S(:, 3) + p(6);
            else
                S(:,:,1) = p(1) * S(:,:,1) + p(4);
                S(:,:,2) = p(2) * S(:,:,2) + p(5);
                S(:,:,3) = p(3) * S(:,:,3) + p(6);
            end;                
        case 'Affine'
            A = reshape(p(1:9), 3, 3);
            if length(size(S)) == 2
                S = S * A' + repmat(p(10:12)', size(S, 1), 1);
            else
                for c = 1:3
                    S(:, :, c) = A(c, 1)*S(:, :, 1)+A(c, 2)*S(:, :, 2)+A(c, 3)*S(:, :, 3)+p(9+c);
                end;
            end;                            
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
end;


