function ng = DIRT_GModel_N(gmodel)
% ng = DIRT_GModel_N(gmodel)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

switch gmodel
    case 'Homography'
        ng = 8;
    case '3PtHomography'
        ng = 6;
    case 'Affine'
        ng = 6;
    case 'Rt'
        ng = 3;
    otherwise
        error(['unknown geometric model ' gmodel]);
end;
