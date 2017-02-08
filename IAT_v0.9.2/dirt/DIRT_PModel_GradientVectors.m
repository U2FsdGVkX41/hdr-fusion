function G = DIRT_PModel_GradientVectors(pmodel, color, Is)
% G = DIRT_PModel_GradientVectors(pmodel, color, Is)
%
% Computes the gradient of the photometric model for the source image.
%
% Inputs:
%  - PModel [string]
%       The photometric model
%  - color [boolean]
%       Indicates color images
%  - Is [column n-vector]
%       The source image at the points of interest
%
% Outputs:
%  - G [n x 2 matrix in the grey-level case and 3n x 2 matrix in the color case]
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~color
    switch pmodel
        case 'GainAndBias'
            G = [ Is ones(length(Is), 1) ];
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
else
    switch pmodel
        case 'SingleGainAndBias'
            G = [ vect(Is') ones(size(Is, 1)*3, 1) ];
        case 'MultipleGainsAndBiases'
            n = size(Is, 1);
            G = zeros(3*n, 6);
            G(1:3:end, 1) = Is(:, 1);
            G(2:3:end, 2) = Is(:, 2);
            G(3:3:end, 3) = Is(:, 3);
            G(1:3:end, 4) = 1;
            G(2:3:end, 5) = 1;
            G(3:3:end, 6) = 1;
        case 'Affine'
            n = size(Is, 1);
            G = zeros(3*n, 12);
            G(1:3:end, 1) = Is(:, 1);
            G(2:3:end, 2) = Is(:, 1);
            G(3:3:end, 3) = Is(:, 1);
            G(1:3:end, 4) = Is(:, 2);
            G(2:3:end, 5) = Is(:, 2);
            G(3:3:end, 6) = Is(:, 2);
            G(1:3:end, 7) = Is(:, 3);
            G(2:3:end, 8) = Is(:, 3);
            G(3:3:end, 9) = Is(:, 3);
            G(1:3:end, 10) = 1;
            G(2:3:end, 11) = 1;
            G(3:3:end, 12) = 1;
        otherwise
            error(['unknown photometric model' pmodel]);
    end;
end;




