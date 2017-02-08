function [warp, mask] = iat_dic(target, source, par)
% [WARP, MASK] = IAT_DIC(IMAGE, TEMPLATE, PAR)
% IAT_DIC is a wrapper of the dual inverse compositional algorithm [1]
% implemented in DIRT package by Adrien Bartoli. 
% (http://isit.u-clermont1.fr/~ab/Research/DIRT_v1p2.zip)
% DIC algorithm extends the inverse compositional scheme 
% by Baker et al. (see the function iat_LucasKanadeIC) to jointly estimate
% a group of photometric transformations along with the warp. 
% IAT_DIC returns only the geometric warp, while the photometric 
% model accounts for any photometric distortions, thus making the
% alignment robust to intensity variations.
% The function relies on pixels around TEMPLATE edges to do the registration.
%
% -->Input:
% IMAGE:                The image that must be finally warped in order to be
%                       similar to TEMPLATE
% TEMPLATE:             The template (fixed) image,
% PAR:                  A struct of parameters with fields:
%
%       PAR.iterations: the number of algorithm iterations (default:100).
%                       The algorithm may terminate earlier owing to a
%                       second criterion used by DIRT
%       PAR.photoModel: the type of photometric transformation.
%                       Valid models:
%                       'GainAndBias', (default for grayscale images)
%                       'SingleGainAndBias' (default for color images),
%                       'MultipleGainsAndBiases',
%                       'Affine'
%       PAR.transform:  the type of geometric transformation.
%                       Valid strings:
%                       'euclidean','affine','homography'
%                       (default: 'homography')
%       PAR.initwarp:   the initial transformation.
%                       Default values:
%                       euclidean, affine: [eye(2) zeros(2,1)]
%                       homography: eye(3)
%       PAR.maskBorder: the size of the TEMPLATE border that is masked 
%                       (shrinkage of ROI) so that warped pixels of IMAGE 
%                       are not mapped outside its borders. 
%                       Default value: 0
%
% -->Output:
% WARP:                 Final estimation of geometric transform
% MASK:                 Binary image with TEMPLATE's pixels used for the
%                       alignemnt.
%
% Notice that the original function DIRT_Registration accepts
% more input parameters, so that a more appropriate tuning is possible for
% each pair of images
% 
% [1] A. Bartoli, Groupwise Geometric and Photometric
% Direct Image Registration, PAMI, vol.30, no.12, 2008
%
% Copyright of DIRT package: Adrien Bartoli, 2008
% For any bugs/questions wrt DIRT, please contact <adrien.bartoli@gmail.com> 
% -------------------
% Author of this wrapper: Georgios Evangelidis,
% Copyright (C) 2014 Georgios Evangelidis
% All rights reserved.
%
% For any bugs/questions wrt the wrapper, please contact <georgios@iatool.net>
% 
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file). It uses functions of
% a third-party package, called DIRT, kindly available by Adrien Bartoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%default parameters
par0.iterations = 100;
par0.transform = 'homography';
if size(target,3)==1
    par0.photoModel = 'GainAndBias';
else
    par0.photoModel = 'SingleGainAndBias';
end
par0.maskBorder = 0;

if exist('par','var')
    if ~isstruct(par)
        error('iat_dic: parameter structure is not a matlab struct');
    end
    if isfield(par,'initwarp') && ~isfield(par,'transform')
        error('iat_dic: when warp is initialized, the type of transform must be defined as well');
    end
    params = iat_merge_param(par0, par);
else
    params = par0;
end

transform = params.transform;
photoModel = params.photoModel;
iter = params.iterations;
maskBorder = params.maskBorder;

if iat_is_transform(transform)
    if (strcmpi(transform,'translation') || ...
            strcmpi(photoModel,'similarity'))
        error('iat_dic: DIC deals only with "Homography", "Affine" and "Euclidean" transforms');
    end
end

if isfield(params,'initwarp')
    ginit = params.initwarp;
    szw = size(ginit);
    switch lower(transform)
        case 'euclidean'
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_dic: warp matrix must be 2x3 for Euclidean transform');
            end
        case 'affine'
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_dic: warp matrix must be 2x3 for affine transform');
            end
        case 'homography'
            if (szw(1)~=3 || szw(2)~=3)
                error('iat_dic: warp matrix must be 3x3 for homography transform');
            end
    end
end

if isfield(params,'initwarp')
    % DIRT uses row/column convention instead of x/y 
    if strcmpi(transform,'homography')
        ginit = [ginit(2,:); ginit(1,:); ginit(3,:)];
        ginit = [ginit(:,2) ginit(:,1) ginit(:,3)];
    else
        ginit = [ginit(2,:); ginit(1,:)];
        ginit = [ginit(:,2) ginit(:,1) ginit(:,3)];
    end
end


if strcmpi(transform,'homography')
    gmodel = 'Homography';
end
if strcmpi(transform,'affine')
    gmodel = 'Affine';
end
if strcmpi(transform,'euclidean')
    gmodel = 'Rt';
end


if exist('photoModel','var')
    if ~(strcmpi(photoModel,'GainAndBias') || ...
            strcmpi(photoModel,'SingleGainAndBias') || ...
            strcmpi(photoModel,'MultipleGainsAndBiases') || ...
            strcmpi(photoModel,'Affine')       )
        error('iat_dic: unknown photometric model');
    end
end

if size(source,3)~=size(target,3)
    error('iat_dic:images should have the same number of channels.');
end

if (size(source,3)==3 && strcmpi(photoModel,'GainAndBias') || ...
        size(source,3)==1 && strcmpi(photoModel,'SingleGainAndBias'))
    error('iat_dic: Use "SingleGainAndBias" with color images and "GainAndBias" with grayscale images');
end


if size(source,3)==3
    sourceGray = rgb2gray(source);
else
    sourceGray = source;
end

% define the region of pixels used for registration.
% This is the area around the edges after removing the borders based on maskBorder
Mask = DIRT_MaskEdges(uint8(sourceGray), 2);
Mask = Mask & DIRT_MaskBorder(maskBorder, size(sourceGray,1), size(sourceGray,2));

mask = Mask;

if ~exist('photoModel','var')
    pre = DIRT_Precomputation(double(source), 'ROI', DIRT_Mask2ROI(Mask), 'gmodel',gmodel, 'no_poptim');
else
    pre = DIRT_Precomputation(double(source), 'ROI', DIRT_Mask2ROI(Mask), 'gmodel',gmodel, 'pmodel', photoModel);
end


% this function does the whole work
if isfield(params,'initwarp')
    reg = DIRT_Registration(double(target), pre, 'max_nb_it', iter, 'ginit', ginit);
else
    reg = DIRT_Registration(double(target), pre, 'max_nb_it', iter);
end

H = reg.g;

if strcmpi(transform,'homography')
    H = H./H(3,3);
    
    % row/color to x/y convention used by IAT
    H = [H(2,:); H(1,:); H(3,:)];
    H = [H(:,2) H(:,1) H(:,3)];
else
    % row/color to x/y convention used by IAT
    H = [H(2,:); H(1,:)];
    H = [H(:,2) H(:,1) H(:,3)];
end

warp = H;

