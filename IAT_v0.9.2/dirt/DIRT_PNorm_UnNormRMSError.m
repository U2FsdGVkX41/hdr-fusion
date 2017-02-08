function e = DIRT_PNorm_UnNormRMSError(params, e)
% e = DIRT_PNorm_UnNormRMSError(params, e)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

e = e * params;