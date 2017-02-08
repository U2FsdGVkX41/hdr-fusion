function [wimage, support] = iat_inverse_warping(image, warp, transform, nx, ny, str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [WIMAGE, SUPPORT] = IAT_INVERSE_WARPING(IMAGE, WARP, TRANSFORM, NX, NY, STR)
% IAT_INVERSE_WARPING implements the inverse warping of IMAGE. The coordinates
% defined by the combinations of NX,NY row-vectors are projected through WARP 
% thus resulting in new subpixel coordinates. The intensities of new 
% coordinates are computed via interpolation of IMAGE.
% For valid interpolation methods see Matlab function INTERP2.
% 
% -->Input:
% IMAGE:                the input image that must be warped,
% WARP:                 the warp transform,
% TRANSFORM:            the type of transformation. Valid string:
%                       {'translation', 'euclidean','similarity','affine','homography'}
% NX:                   the x-coordinate values of horizontal side of ROI
%                       (i.e. [xmin:xmax]),
% NY:                   the y-coordinate values of vertical side of ROI
%                       (i.e. [ymin:ymax]),
% STR:                  (optional input) the string that corresponds to 
%                       interpolation method:
%                       'linear', 'cubic' etc (for details see the
%                       Matlab function INTERP2), Default: 'linear'
%
% 
% -->Output:
% WIMAGE:               the warped (interpolated) image
% SUPPORT:              a binary double image that that shows which areas 
%                       in WIMAGE come from the support area of 
%                       IMAGE and which from outside the borders
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

if nargin==5
    str = 'linear';
end

if ~iat_is_transform(transform)
    error('iat_inverse_warping: non-valid transform name');
end


% in affine, similarity and euclidean case, make the warp matrix 3x3
if (strcmpi(transform,'affine')||strcmpi(transform,'euclidean') || strcmpi(transform,'similarity'))
    if size(warp,1)==2
        warp=[warp; 0 0 1];
    end
end

% in translation case, make the warp matrix 3x3
if strcmp(transform,'translation')
    if size(warp,2)==1 && size(warp,1)==2
        warp = [eye(2) warp];
        warp = [warp; 0 0 1];
    end
    if size(warp,2)==3 && size(warp,1)==2
        warp = [warp; 0 0 1];
    end
    
end

[xx,yy] = meshgrid(nx, ny);
xy=[xx(:)';yy(:)';ones(1,length(yy(:)))];

%3x3 matrix transformation
A = warp;
if strcmp(transform,'homography')
    A = A./A(3,3);
end

% new coordinates
xy_prime = A * xy;

if strcmp(transform,'homography')
    
    % division due to homogeneous coordinates
    xy_prime(1,:) = xy_prime(1,:)./xy_prime(3,:);
    xy_prime(2,:) = xy_prime(2,:)./xy_prime(3,:);
end

% Ignore third row
xy_prime = xy_prime(1:2,:);

image = double(image);
wimage = zeros(length(ny),length(nx),size(image,3));
for i = 1:size(image,3)
    
    % Subpixel interpolation
    tempimage = interp2(image(:,:,i), xy_prime(1,:), xy_prime(2,:), str);
    if i==1
        support = double(~isnan(tempimage));
        support = reshape(support,length(ny),length(nx));
    end
    tempimage(isnan(tempimage))=0;%replace Nan
    wimage(:,:,i) = reshape(tempimage,length(ny),length(nx));
    
end