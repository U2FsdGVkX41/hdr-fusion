function [siftImage, gridX, gridY] = iat_dense_sift(im, patch_size, grid_spacing, varargin)
% [SIFTIMAGE, GRIDX,GRIDY] = IAT_DENSE_SIFT(IMAGE, PS, GSPACING)
% IAT_DENSE_SIFT creates SIFTIMAGE from IMAGE, when the pixels of latter
% are replaced by SIFT descriptors [1] that describe their surrounding
% area od size PSxPS. SIFTIMAGE is defined on the meshgrid [GRIDX,GRIDy]
% which results from IMAGE's sampling with GRIDSPACING factor (and
% after ignoring PS/2 rows and columns from all image sides).
%
% Input arguments:
% IMAGE:                the input image
% PS:                   The width (=height) of the square patch
%                       that is described by SIFT vectors.
% GRIDSPACING:          the donwscale factor that defines the resolution of
%                       SIFTIMAGE
%
%
% -->Optional Paramaters:
% [SIFTIMAGE, GRIDX,GRIDY] = iat_dense_sift(IMAGE, PS, GSPACING, 'PARAM1', PARAM1VALUE,...)
% The user can define his own parameter values insted of default ones.
% These parameters are:
%
% 'numAngles':          The quantization step (angle) for gradient orientation
%                       (default: 8)
%
% 'numBins':            The number of spatial bins in each dimension. The
%                       size of each SIFT descriptor is numAngles*numBins*numBins
%                       (default: 4)
%
% 'alpha':              The attenuation of angles; it must be odd number
%                       (default: 9)
%
% 'sigma':              The scale of gaussian kernel for computing
%                       DOG (default: 1). When sigma is scalar, the size of
%                       kernel is (4*ceil(sigma)+1)X(4*ceil(sigma)+1). When
%                       sigma is a 2-element vector, i.e. [sigmaX, sigmaY],
%                       the size of the kernel is
%                       (4*ceil(sigmaY)+1)X(4*ceil(sigmaX)+1)
%
% Output arguments
% SIFTIMAGE:            the output image
% GRIDX & GRIDY:        the meshgrid of SIFTIMAGE's support area.
%
%
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

im = double(im);
im = mean(im,3);
im = im /max(im(:));

% input parser
par = inputParser;
addOptional(par,'numAngles',8);
addOptional(par,'numBins',4);
addOptional(par,'alpha',9);
addOptional(par,'sigma',1);

% merge default with provided values
parse(par,varargin{:});

% get the parameters
num_angles = par.Results.numAngles;
num_bins = par.Results.numBins;
alpha = par.Results.alpha;
sigma_edge = par.Results.sigma;

if (mod(alpha,2)~=1)
    error('iat_dense_sift: Parameter "alpha" must be odd number');
end

angle_step = 2 * pi / num_angles;
angles = 0:angle_step:2*pi;
angles(num_angles+1) = []; % bin centers

[height, width] = size(im);


[GX,GY] = gradient(gen_gauss(sigma_edge));
GX = GX * 2 ./ sum(sum(abs(GX)));
GY = GY * 2 ./ sum(sum(abs(GY)));
I_X = filter2(GX, im, 'same'); % vertical edges
I_Y = filter2(GY, im, 'same'); % horizontal edges
I_mag = sqrt(I_X.^2 + I_Y.^2); % gradient magnitude
I_theta = atan2(I_Y,I_X);
I_theta(isnan(I_theta)) = 0;

% grid of SIFTIMAGE
gridX = patch_size/2+1:grid_spacing:width-patch_size/2;
gridY = patch_size/2+1:grid_spacing:height-patch_size/2;

% orientation images
I_ori = zeros([height, width, num_angles], 'single');

% for each histogram angle
cosI = cos(I_theta);
sinI = sin(I_theta);
for a=1:num_angles
    % compute each orientation channel
    tmp = (cosI*cos(angles(a))+sinI*sin(angles(a))).^alpha;
    tmp = tmp .* (tmp > 0);
    
    % weight by magnitude
    I_ori(:,:,a) = tmp .* I_mag;
end

% Convolution formulation:
r = patch_size/2;
cx = r - 0.5;
sample_res = patch_size/num_bins;
weight_x = abs((1:patch_size) - cx)/sample_res;
weight_x = (1 - weight_x) .* (weight_x <= 1);


for a = 1:num_angles
    I_ori(:,:,a) = conv2(weight_x, weight_x', I_ori(:,:,a), 'same');
end

% Sample SIFT bins at valid locations (without boundary artifacts)
% find coordinates of sample points (bin centers)
[sample_x, sample_y] = meshgrid(linspace(1,patch_size+1,num_bins+1));
sample_x = sample_x(1:num_bins,1:num_bins);
sample_x = sample_x(:)-patch_size/2;
sample_y = sample_y(1:num_bins,1:num_bins);
sample_y = sample_y(:)-patch_size/2;


siftImage = zeros([length(gridY) length(gridX) num_angles*num_bins*num_bins], 'single');
b = 0;
for n = 1:num_bins*num_bins
    siftImage(:,:,b+1:b+num_angles) = I_ori(round(gridY+sample_y(n)), round(gridX+sample_x(n)), :);
    b = b+num_angles;
end

% Outputs:
[gridX,gridY] = meshgrid(gridX, gridY);
[nrows, ncols, cols] = size(siftImage);

% normalize SIFT descriptors
siftImage = reshape(siftImage, [nrows*ncols num_angles*num_bins*num_bins]);
siftImage = normalize_sift(siftImage);
siftImage = reshape(siftImage, [nrows ncols num_angles*num_bins*num_bins]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = pad(x, D)

[nrows, ncols, cols] = size(siftImage);
hgt = nrows+2*D;
wid = ncols+2*D;
PADVAL = 0;

x = [repmat(PADVAL, [hgt Dx cols]) ...
    [repmat(PADVAL, [Dy ncols cols]); x; repmat(PADVAL, [Dy-1 ncols cols])] ...
    repmat(PADVAL, [hgt Dx-1 cols])];


function G=gen_gauss(sigma)
% generate a gauss kernel

if all(size(sigma)==[1, 1])
    % isotropic gaussian
    f_wid = 4 * ceil(sigma) + 1;
    G = fspecial('gaussian', f_wid, sigma);
    
else
    % anisotropic gaussian
    f_wid_x = 2 * ceil(sigma(1));
    f_wid_y = 2 * ceil(sigma(2));
    G_x = normpdf(-f_wid_x:f_wid_x,0,sigma(1));
    G_y = normpdf(-f_wid_y:f_wid_y,0,sigma(2));
    G = G_y' * G_x;
    G = G/sum(G(:));
end


function sift_arr = normalize_sift(sift_arr)
% normalize SIFT descriptors of input sift_arr (each row is a descriptor)

% find indices of descriptors to be normalized (those whose norm is larger than 1)
tmp = sqrt(sum(sift_arr.^2, 2));
normalize_ind = (tmp > 1);

sift_arr_norm = sift_arr(normalize_ind,:);
sift_arr_norm = sift_arr_norm ./ repmat(tmp(normalize_ind,:), [1 size(sift_arr,2)]);

% suppress large gradients
sift_arr_norm(sift_arr_norm > 0.2) = 0.2;

% finally, renormalize to unit length
tmp = sqrt(sum(sift_arr_norm.^2, 2));
sift_arr_norm = sift_arr_norm ./ repmat(tmp, [1 size(sift_arr,2)]);

sift_arr(normalize_ind,:) = sift_arr_norm;