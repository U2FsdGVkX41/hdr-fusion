function DIRT_GModel_Disp(gmodel, g)
% DIRT_GModel_Disp(gmodel, g)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

switch gmodel
    case {'Homography','3PtHomography','Affine','Rt'}
        fprintf('');
    otherwise
        error(['unknown geometric model ' gmodel]);
end;
