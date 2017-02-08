function [imdiff, grayerror] = iat_error2gray( im1, im2, support)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [IMDIFF, GRAYERROR] = IAT_ERROR2GRAY(IM1, IM2, SUPPORT)
% IAT_ERROR2GRAY returns the difference IM1-IM2 in IMDIFF and its 
% translation to a grayscale image in GRAYERROR. Both take into account 
% the mask SUPPORT when provided, i.e. error is computed in the area 
% marked by binary image SUPPORT. Color images are converted to 
% grayscale before the subtraction
%
% -->Input:
% IM1:                  The first image
% IM2:                  The second image. The images are scaled so that the 
%                       intensities are in range [0,255] before the
%                       subtraction
% SUPPORT:              The support area for which the error is computed
%                       (default: all-one mask of size(IM1))
%
% -->Output:
% IMDIFF:               The image difference.
% GRAYERROR:            The image difference as UINT8 (grascale) image. 
%
% -------------------
% Authors: Georgios Evangelidis, Panagiotis Anatolitis
% Copyright (C) 2013 Georgios Evangelidis
% All rights reserved.
%
% For any bugs, please contact <georgios.evangelidis@inria.fr> or
% <anatolitis@ceid.upatras.gr>
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==2
    support = ones(size(im1,1), size(im1,2));
end

if (size(im1,3)~=size(im2,3))
    error('iat_error2gray: images may not have the same color format');
end

if ~all(size(im1)==size(im2))
    error('iat_error2gray: images must have the same size');
end

im1 = double(im1);
im2 = double(im2);

[~,~,K]=size(im1);

[sr, sc] = size(support);

if (size(im1,1)~=sr) || (size(im1,2)~=sc)
    error('iat_error2gray: inconsistent size between images and support mask');
end


if max(im1(:))<(1+eps)
    im1 = im1*255;
end

if max(im2(:))<(1+eps)
    im2 = im2*255;
end

if K==3
    im1=double(rgb2gray(uint8(im1)));
    im2=double(rgb2gray(uint8(im2)));
end

im1 = im1.* support;
im2 = im2.* support;

imdiff = im1-im2;

grayerror = abs(imdiff);
 
% grayerror = grayerror-min(grayerror(:));
% grayerror = grayerror/max(grayerror(:));
% grayerror = grayerror*255;

grayerror = uint8(255-grayerror);

end

