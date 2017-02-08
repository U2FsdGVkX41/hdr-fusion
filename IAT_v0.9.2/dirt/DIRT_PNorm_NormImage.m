function I = DIRT_PNorm_NormImage(params, I)
% I = DIRT_PNorm_NormImage(params, I)
%
% Note: the normalization must be such that it can be applied after image
% derivation.
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

I = I / params;
