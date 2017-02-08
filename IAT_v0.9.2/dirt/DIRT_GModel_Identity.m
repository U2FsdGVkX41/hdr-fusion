function g = DIRT_GModel_Identity(gmodel,gopt)
% g = DIRT_GModel_Identity(gmodel,gopt)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

switch gmodel
    case 'Homography'
        g = eye(3);
    case '3PtHomography'
        g = [ eye(2) gopt.tq-gopt.sq ; 0 0 1 ];
    case {'Affine','Rt'}
        g = [ eye(2) zeros(2,1) ];
    otherwise
        error(['unknown geometric model ' gmodel]);
end;
