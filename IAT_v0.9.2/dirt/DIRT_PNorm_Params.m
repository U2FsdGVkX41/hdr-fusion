function params = DIRT_PNorm_Params(flag, sI)
% params = DIRT_PNorm_Params(flag, sI)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~flag
    params = 1;
else
    params = max(max(max(sI)));
end;