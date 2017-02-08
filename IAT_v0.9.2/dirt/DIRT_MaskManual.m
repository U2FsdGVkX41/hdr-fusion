function [M, xi, yi] = DIRT_MaskManual(sI)
% [M, xi, yi] = DIRT_MaskManual(sI)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

[M, xi, yi] = roipoly(uint8(sI));
if nargout > 1 oxi = xi; end;
if nargout > 2 oyi = yi; end;
close;