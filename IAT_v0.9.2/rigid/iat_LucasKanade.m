function [warp] = iat_LucasKanade(image, template, par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [WARP] = iat_LucasKanade(IMAGE, TEMPLATE, PAR)
% iat_LucasKanade implements the forwards-addiitve version of Lucas-Kanade
% image alignment algorithm [1]. It computes the geometric transformation 
% that should be applied to IMAGE in order to get an image similar to 
% TEMPLATE. The optimum transformation is the one that minimizes the
% squared image difference between TEMPLATE and warped IMAGE.
%
%
% -->Input:
% IMAGE:                The image that must be warped in order to be
%                       similar to TEMPLATE
% TEMPLATE:             The target image,
% PAR:                  A struct of parameters with fields:
%
%       PAR.iterations: the number of algorithm's iteration (defualt:50)
%       PAR.levels:     the number of levels for multi-resolution axecution
%                       (default: 1)
%       PAR.transform:  the type of geometric transformation. Valid strings:
%                       'translation','euclidean','affine','homography'
%                       (default: 'affine')
%       PAR.initwarp:   the initial transformation. Default values:
%                       translation: zeros(2,1)
%                       euclidean: [eye(2) zeros(2,1)]
%                       affine: [eye(2) zeros(2,1)]
%                       homography: eye(3)
%
%
% -->Output:
% WARP:                 The final estimated transformation
%
% References:
% [1] S. Baker, I. Matthews, "Lucas-Kande 20 years on: A unifying framework,
% Part I", IJCV, vol.56, no. 3, 2004
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<2
    error('iat_LucasKanade:Not enough input arguments');
end

%default parameters
par0.iterations = 50;
par0.levels = 1;
par0.transform = 'affine';

if exist('par','var')
    if ~isstruct(par)
        error('iat_LucasKanade: the datatype of parameters is not a matlab struct');
    end
    
    if isfield(par,'initwarp') && ~isfield(par,'transform')
        error('iat_LucasKanade: when you initialize the warp, you should define the type of transform as well');
    end
    params = iat_merge_param(par0, par);
else
    params = par0;
end

if ~iat_is_transform(params.transform)
    error('iat_LucasKanade: unknown transform type. Check the field .transform in parameters structure');
end

if strcmpi(params.transform,'similarity')
    params.transform = 'affine';
    warning('iat_LucasKanade: Lukas-Kanade does not support similarity transform. Warp was automatically changed to affine')
end
    
transform = params.transform;    



if isfield(params,'initwarp')
    warp = params.initwarp;
    szw = size(warp);
    switch lower(transform)
        case 'translation'
            nop =2;
            if (szw(1)~=2 || szw(2)~=1)
                error('iat_LucasKanade: the warp matrix must be 2x1 for translation transform');
            end
        case 'euclidean'
            nop = 3;
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_LucasKanade: the warp matrix must be 2x3 for euclidean transform');
            end
        case 'affine'
            nop=6;
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_LucasKanade: the warp matrix must be 2x3 for affine transform');
            end
        case 'homography'
            nop = 8;
            if (szw(1)~=3 || szw(2)~=3)
                error('iat_LucasKanade: the warp matrix must be 3x3 for homography transform');
            end
    end
else
    switch lower(transform)
        case 'translation'
            warp = zeros(2,1);
            nop =2;
        case 'euclidean'
            nop = 3;
            warp = [eye(2) zeros(2,1)];
        case 'affine'
            nop=6;
            warp = [eye(2) zeros(2,1)];
        case 'homography'
            nop = 8;
            warp = eye(3);
    end
    
end


levels = params.levels;
noi = params.iterations;

% Number of color channels for both image and template
sZi3 = size(image,3);
sZt3 = size(template,3);

% Color format validity check for image (RGB or gray-scale)
if sZi3>1
    if ((sZi3==2) || (sZi3>3))
        error('iat_LucasKanade: Unknown color image format: check the number of channels');
    else
        image=rgb2gray(uint8(image));
    end
end

% Color format validity check for image (RGB or gray-scale)
if sZt3>1
    if ((sZt3==2) || (sZt3>3))
        error('iat_LucasKanade: Unknown color image format: check the number of channels');
    else
        template = rgb2gray(uint8(template));
    end
end

% Converting template and image to doubles
template = double(template);
image = double(image);

%% pyramid images
% The following for-loop creates pyramid images in cells IM and TEMP with varying names
% The variables IM{1} and TEMP{1} are the images with the highest resoltuion

TEMP{1} = template;
IM{1} = image;

% Enable smoothing (optional)
% f = fspecial('gaussian',[7 7],.5);
% TEMP{1} = imfilter(template,f);
% IM{1} = imfilter(image,f);

for nol=2:levels
    IM{nol} = imresize(IM{nol-1},.5);
    TEMP{nol} = imresize(TEMP{nol-1},.5);
end

% in case of pyramid implementation, the initial transformation must be
% appropriately modified
for ii=1:levels-1
    warp=iat_warp_updown(warp, transform, 0);
end



if levels==1
    disp('Lucas-Kanade is running in single-level mode....');
else
    disp('Lucas-Kanade is running in multi-level mode....');
end

%% Run Lucas-Kanade algorithm for each level of pyramid
for nol=levels:-1:1
        if levels>1
        fprintf('Level %d...', nol);
    end
    im = IM{nol};
    [vx,vy]=gradient(im);
    
    temp = TEMP{nol};
    
   [A,B]=size(temp);
    margin = 0; % no margin (enable a margin if you want)
    
    nx=margin+1:B-margin;
    ny=margin+1:A-margin;
    temp=double(temp(ny,nx,:));
    
    
    for i=1:noi
        
        %disp(['LucasKanade: Level: ' num2str(nol) ', Iteration: ' num2str(i)])
        %Image interpolation method
        str='linear'; % bilinear interpolation 
        %str='cubic'; % cubic ibterpolation
        
        wim = iat_inverse_warping(im, warp, transform, nx, ny, str); %inverse (backward) warping
         
        if (i == noi) % the algorithm is executed (noi-1) times
            break;
        end
        
        % Gradient Image interpolation (warped gradients)
        wvx = iat_inverse_warping(vx, warp, transform, nx, ny, str);
        wvy = iat_inverse_warping(vy, warp, transform, nx, ny, str);
        
        % Compute the jacobian of warp transform
        J = iat_warp_jacobian(nx, ny, warp, transform);
        
        % Compute the jacobian of warped image wrt parameters (steepest
        % descent image)
        G = iat_image_jacobian(wvx, wvy, J, nop);
        
        % Compute Hessian and its inverse
        C= G' * G;% C: Hessian matrix
        %i_C = inv(C);
               
        % Compute error vector
        imerror = temp - wim;
        
        % Compute the projection of error vector into Jacobian G
        Ge = G' * imerror(:);
        
        % Compute the optimum parameter correction vector
        delta_p = C\Ge;
        
        % Update parmaters
        warp = iat_warp_update(warp, delta_p, transform);
        
        
    end
    
    % modify the parameteres appropriately for next pyramid level
    if (nol>1)
        warp = iat_warp_updown(warp, transform,1);
    end
      if levels>1
    fprintf('Done\n');
    end  
end
