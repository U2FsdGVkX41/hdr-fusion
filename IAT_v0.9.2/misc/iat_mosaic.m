function [ mosaic ] = iat_mosaic( img1, img2, H )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ MOSAIC ] = IAT_MOSAIC( IMAGE1, IMAGE2, H )
% IAT_MOSAIC creates a mosaic image given two input images and the
% transformation matrix that relates them. IMAGE1 is the reference image
% and is stiched with IMAGE2 after the transformation of the latter by H. 
% H should be given as a 3x3 array. For example, give [1 0 tx; 0 1 ty; 0 0 1] 
% for 2D translation.
%
% -->Input:
% IMAGE1:               The first image
% IMAGE2:               The second image
% H:                    The 3x3 transformation matrix so that 
%                       IMAGE2(H(x,y))->IMAGE1(x,y)
%
% -->Output:
% MOSAIC:               The final (mosaic) image.
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

img1=double(img1);
img2=double(img2);

[xA, yA, zA]=size(img1);
[xB, yB, zB]=size(img2);

if zA~=zB
    error('iat_mosaic: images must have the same number of channels');
end


if size(H,1)==2 && size(H,2)==3 % affine, similarity, euclidean
    H = [H; 0 0 1];
end

[m,n]=size(H);

if m~=3 || n~=3
    error('iat_mosaic: Transform should be given as 3x3 or 2x3 array (e.g. [1 0 tx; 0 1 ty; 0 0 1] for 2D translation)');
end

% normalize homography elements with H(3,3)
H=H./H(3,3);


% Send corners of B to A's frame
cornersB2A = H\[1 1 yB yB;1 xB 1 xB;1 1 1 1];

% Get non-homogeneous coordinates
cornersB2A = iat_remove_scale(cornersB2A);

% estimate the grid to be warped
min_x = round(min([cornersB2A(1,:) 1]));
max_x = round(max([cornersB2A(1,:) yA]));

min_y = round(min([cornersB2A(2,:) 1]));
max_y = round(max([cornersB2A(2,:) xA]));


% inverse-warp of image B
[mosaic, support] = iat_inverse_warping( img2, H, 'homography', min_x:max_x, min_y:max_y );

mosaic(~support) = 255;

if size(mosaic,3)==3
    tempm = mosaic(:,:,2);
    tempm(~support) = 255;
    mosaic(:,:,2) = tempm;
    tempm = mosaic(:,:,3);
    tempm(~support) = 255;
    mosaic(:,:,3) = tempm;
end
    

%position of A's corners in the mosaic
A_X_min = abs(min_x-1)+1;
A_X_max = A_X_min+yA-1;

A_Y_min = abs(min_y-1)+1;
A_Y_max = A_Y_min+xA-1;

% superimpose image A
mosaic(A_Y_min:A_Y_max, A_X_min:A_X_max,:) = img1;


end

