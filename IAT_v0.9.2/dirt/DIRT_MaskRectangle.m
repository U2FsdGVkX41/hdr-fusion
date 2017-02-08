function M = DIRT_MaskRectangle(sI, r, c, rHalfSize, cHalfSize)
% M = DIRT_MaskRectangle(sI, r, c, rHalfSize, cHalfSize)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

M = zeros(size(sI, 1), size(sI, 2));
M(r-rHalfSize:r+rHalfSize, c-cHalfSize:c+cHalfSize) = 1;
