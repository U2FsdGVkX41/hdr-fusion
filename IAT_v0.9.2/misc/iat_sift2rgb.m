function img=iat_sift2rgb(sift_img)
% IM = IAT_SIFT2RGB(SIFTIM)
% IAT_SIFT2RGB maps an 128-channel SIFT-image SIFTIMG to RGB (3-channel) 
% space using a pretrained PCA model, so that a RGB image IMG is obtained. 
% After projecting SIFT pixels to the three most principal components, 
% the RGB channels are obtained by inverting the following mapping:
% 
% first component  -> R+G+B
% second component -> R-G
% third component  -> (R+G)/2-B
% 
% -->Input
% SIFTIM:               The input SIFT-image
%
% -->Output
% IM:                   The color (RGB) image
%
% -------------------
% Authors: Ce Liu, Georgios Evangelidis
% Copyright (C) 2013 Ce Liu
% All rights reserved.
%
% This function is a minor modification of the original function written by
% Ce Liu. For any bugs, please contact <celiu@microsoft.com> or
% <georgios@iatool.net>
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[nrows, ncols, nfeatures] = size(sift_img);

if nfeatures~=128
    error('iat_sift2rgb: SIFT-image must be of size MxNx128');
end

load ('pcSIFT', 'pcSIFT')

SIFTpca = pcSIFT(:,1:3)'*double(reshape(sift_img, [nrows*ncols nfeatures]))';

A = inv([1 1 1; 1 -1 0; .5 .5 -1]);

img = A * SIFTpca(1:3,:);
img = reshape(img', [nrows ncols 3]);
img = img - min(img(:));
img = img / max(img(:));
img = uint8(255*img);
