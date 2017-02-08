function g = DIRT_GModel_Init(gmodel,sR,sC,tR,tC)
% g = DIRT_GModel_Init(gmodel,sR,sC,tR,tC)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

switch gmodel
    case {'Homography','3PtHomography'}
        g = [
            tR/sR 0 0
            0 tC/sC 0
            0 0 1
            ];
    case 'Affine'
        g = [
            tR/sR 0 0
            0 tC/sC 0
            ];
    case 'Rt'
        error('cannot initialize the Rt model from the image sizes');
    otherwise
        error(['unknown geometric model ' gmodel]);
end;
