function M = DIRT_MaskEdges(sI, r)
% M = DIRT_MaskEdges(sI, r)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~exist('r', 'var') r = 2; end;

if ABIsRGB(sI)
    error('A grey-level image is required');
end;
M = edge(sI, 'prewitt');
M = imdilate(M, strel('disk', r));
