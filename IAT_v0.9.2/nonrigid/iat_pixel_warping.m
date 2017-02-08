function [ wimg, support ] = iat_pixel_warping( img, dx, dy )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [WIMG, SUPPRORT] =  = IAT_PIXEL_WARPING(IMG,DX,DY)
%
% IAT_PIXEL_WARPING performs an inverse warping in a pixel basis, i.e., 
% every pixel is warped separately according to horizontal and vertical 
% displacements DX and DY respectively.
%
% -->Input:
% IMG:                  The input image
% DX:                   the array with horizontal displacements.
% DY:                   the array with vertical displacements. Pixel (i,j)
%                       is horisontally translated by DX(i.j) and vertically 
%                       by DY(i,j)
%
% -->Output:
% WIMG:                 The warped image
% SUPPORT:              a binary image that shows which areas 
%                       of WIMG come from the support area of 
%                       IMAGE and which from outside the borders
%
% -------------------
% Authors: Ce Liu, Georgios Evangelidis
% Copyright (C) 2013 Ce Liu
% All rights reserved.
%
% This function is a modification of the original function written by
% Ce Liu. For any bugs, please contact <celiu@microsoft.com> or
% <georgios@iatool.net>
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input Checking
[rim, cim, zim]=size(img);
[rx, cx]=size(dx);
[ry, cy]=size(dy);

if (rim~=rx) || ( rim~=ry) || (cim~=cx) || (cim~=cy)
    error('iat_pixel_warping: Image and displacement arrays must have the same size');
end

% Create image and velocity grid
[xx,yy]=meshgrid(1:cim,1:rim);
XX = xx;
YY = yy;

% Adding displacements to pixels
XX=XX+dx;
YY=YY+dy;

% Mask that shows which pixels stay in the support area
support=XX<1 | XX>cim | YY<1 | YY>rim;

wimg = zeros(size(img));
img  = double(img);

% Interpolate pixel values to the new positions
for i=1:zim
    temp=interp2(xx,yy,img(:,:,i),XX,YY,'bicubic');
    temp(support)=0;
    wimg(:,:,i)= temp;
end

support = 1-support;           

if isa(img,'uint8')
    wimg = uint8(wimg);
end
