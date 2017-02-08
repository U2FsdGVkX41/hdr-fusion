function ROI = DIRT_Mask2ROI(M)
% ROI = DIRT_Mask2ROI(M)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

ROI.Rind = find(M);
ROI.Dind = find(~M);
ROI.Rn   = length(ROI.Rind);
