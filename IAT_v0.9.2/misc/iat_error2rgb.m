function [imerr] = iat_error2rgb( im1, im2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMERR = IAT_ERROR2RGB(IM1, IM2)
% IAT_ERROR2RGB returns an RGB image that mixes the channels of IM1 and IM2.
% When input images are grayscale, then R,B, channels come from IM1 and 
% G channel comes from IM2. When it coms to grascale images, R,B channel 
% are equal to IM1 and G channel is equal to IM2. This way, misalignments 
% between IM1 and IM2 appear with hot pink and lawn green colors.
%
% -->Input:
% IM1:                  The first image
% IM2:                  The second image (of the same size with IM1)
%
% -->Output:
% IMERR:                the RGB image that visualizes misalignments between
%                       IM1 and IM2
%
% -------------------
% Authors: Georgios Evangelidis
% Copyright (C) 2013 Georgios Evangelidis
% All rights reserved.
%
% For any bugs, please contact <georgios@iatool.net> or
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (size(im1,3)~=size(im2,3))
    error('iat_error2rgb: images must have the same number of channels');
end

if ~all(size(im1)==size(im2))
    error('iat_error2rgb: images must have the same size');
end

if ~strcmp(class(im1),class(im2))
    error('iat_error2rgb: images must be objects of the same class (e.g. uint8)');
end


[row,col,K]=size(im1);

imerr = zeros(row,col,3);

if K==1
    imerr = cat(3,im1,im2,im1);
elseif K==3
    imerr = cat(3,im1(:,:,1),im2(:,:,2),im1(:,:,3));
else
    error('iat_error2rgb: unknown color format');
end

if max(imerr(:))<=1
    imerr = imerr*255;
end

imerr = uint8(imerr);

