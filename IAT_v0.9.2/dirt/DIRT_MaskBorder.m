function M = DIRT_MaskBorder(bsz, nr, nc)
% M = DIRT_MaskBorder(bsz, nr, nc)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

M = zeros(nr, nc);
M(bsz+1:nr-bsz, bsz+1:nc-bsz) = 1;

