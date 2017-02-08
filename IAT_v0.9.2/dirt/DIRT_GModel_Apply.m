function [ tR, tC ] = DIRT_GModel_Apply(gmodel, g, sR, sC)
% [ tR, tC ] = DIRT_GModel_Apply(gmodel, g, sR, sC)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

switch gmodel
    case {'Homography','3PtHomography'}
        t3 = sR*g(3, 1)+ sC*g(3, 2) + g(3, 3);
        tR = (sR*g(1, 1) + sC*g(1, 2) + g(1, 3)) ./ t3;
        tC = (sR*g(2, 1) + sC*g(2, 2) + g(2, 3)) ./ t3;
    case {'Affine','Rt'}
        tR = sR*g(1, 1) + sC*g(1, 2) + g(1, 3);
        tC = sR*g(2, 1) + sC*g(2, 2) + g(2, 3);
    otherwise
        error(['unknown geometric model ' gmodel]);
end;
