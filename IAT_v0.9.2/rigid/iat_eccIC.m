function [warp, rho] = iat_eccIC(image, template, par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [WARP, RHO] = iat_eccIC(IMAGE, TEMPLATE, PAR)
% iat_eccIC implements the inverse-compositional version of ECC image alignment 
% algorithm [1]. It computes the geometric transformation that should be 
% applied to IMAGE in order to get an image similar to TEMPLATE. The 
% optimum transformation is the one that maximizes the Enhanced Correlation
% Coefficient [1] between TEMPLATE and warped IMAGE. It is more efficient
% than the forwards-additive version (function: iat_ecc) at the cost of
% lower accuracy (see [1] for comparison)
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
% RHO:                  The final correlation coefficient between template
%                       and final warped image
%
% Reference: 
% [1] "G.D.Evangelidis, E.Z.Psarakis, Parametric Image Alignment using
% Enhanced Correlation Coefficient, IEEE Trans. on PAMI, vol.30, no.10, 2008"
%
% -------------------
% Authors: Georgios Evangelidis
% Copyright (C) 2013 Georgios Evangelidis
% All rights reserved.
%
% For any bugs, please contact <georgios.evangelidis@inria.fr> or
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<2
    error('iat_ecc_ic:Not enough input arguments');
end

%default parameters
par0.iterations = 50;
par0.levels = 1;
par0.transform = 'affine';

if exist('par','var')
    if ~isstruct(par)
        error('iat_ecc_ic: the datatype of parameters is not a matlab struct');
    end
    
    if isfield(par,'initwarp') && ~isfield(par,'transform')
        error('iat_ecc_ic: when you initialize the warp, you should define the type of transform as well');
    end
    params = iat_merge_param(par0, par);
else
    params = par0;
end

if ~iat_is_transform(params.transform)
    error('iat_ecc_ic: unknown transform type. Check the field .transform in parameters structure');
end

if strcmpi(params.transform,'similarity')
    params.transform = 'affine';
    warning('iat_ecc_i: ECC-IC does not support similarity transform. Warp was automatically changed to affine')
end
    
transform = params.transform;    
    
if isfield(params,'initwarp')
    warp = params.initwarp;
    szw = size(warp);
    switch lower(transform)
        case 'translation'
            nop =2;
            if (szw(1)~=2 || szw(2)~=1)
                error('iat_ecc_ic: the warp matrix must be 2x1 for translation transform');
            end
        case 'euclidean'
            nop = 3;
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_ecc_ic: the warp matrix must be 2x3 for euclidean transform');
            end
        case 'affine'
            nop=6;
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_ecc_ic: the warp matrix must be 2x3 for affine transform');
            end
        case 'homography'
            nop = 8;
            if (szw(1)~=3 || szw(2)~=3)
                error('iat_ecc_ic: the warp matrix must be 3x3 for homography transform');
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


break_flag=0;

levels = params.levels;
noi = params.iterations;

% Number of color channels for both image and template
sZi3 = size(image,3);
sZt3 = size(template,3);


% Color format validity check for image (RGB or gray-scale)
if sZi3>1
    if ((sZi3==2) || (sZi3>3))
        error('iat_ecc_ic: Unknown color image format: check the number of channels');
    else
        image=rgb2gray(uint8(image));
    end
end

% Color format validity check for image (RGB or gray-scale)
if sZt3>1
    if ((sZt3==2) || (sZt3>3))
        error('iat_ecc_ic: Unknown color image format: check the number of channels');
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
    disp('ECC-IC is running in single-level mode....');
else
    disp('ECC-IC is running in multi-level mode....');
end

%% Run ECC algorithm for each level of pyramid
for nol=levels:-1:1
    if levels>1
        fprintf('Level %d...', nol);
    end
    
    im = IM{nol};
    
    temp = TEMP{nol};
    
    [vx,vy]=gradient(temp);
     
    [A,B]=size(temp);
    % Warning for tiny images
    if prod([A,B])<400
        disp(' -> ECC-IC Warning: The size of images in high pyramid levels is quite small and it may cause errors.');
        disp(' -> Try fewer levels or larger images to avoid such errors.');
        disp(' -> Press any key to continue.')
        pause
    end
    
    margin = 0; % no margin
    
    % Uncomment the follwing lines if you want to consider margins
    % m0 = mean([A,B]);
    % margin = floor(m0*.05/(2^(nol-1))); % 5-percent of the mean of [height,width]
    
    nx=margin+1:B-margin;
    ny=margin+1:A-margin;
    temp=double(temp(ny,nx,:));
    
    
    % Pre-computations
    
    % Compute the jacobian of identity transform
    J = iat_warp_jacobian(nx, ny, eye(3), transform);
    
    % Compute the jacobian of warped image wrt parameters (matrix G in the paper)
    G = iat_image_jacobian(vx, vy, J, nop);
    
    % Compute Hessian and its inverse
    C= G' * G;% C: Hessian matrix
    
    
    %% ECC, Forwards Additive Algorithm -------------------------------
    for i=1:noi
        
        %disp(['ECC: Level: ' num2str(nol) ', Iteration: ' num2str(i)])
        %Image interpolation method
        str='linear'; % bilinear interpolation
        %str='cubic'; % cubi ibterpolation
        
        [wim, ones_map] = iat_inverse_warping(im, warp, transform, nx, ny, str); %inverse (backward) warping
        
        % consider the overap for the zero-mean process
        numOfElem = sum(ones_map(:)~=0);
        meanOfWim = sum(wim(ones_map~=0))/numOfElem;
        meanOfTemp = sum(temp(ones_map~=0))/numOfElem;
        
        % Compute zero-mean images
        wim = wim-meanOfWim;% zero-mean image; is useful for brightness change compensation, otherwise you can comment this line
        tempzm = temp-meanOfTemp; % zero-mean template
        
        wim(ones_map==0) = 0; % reject pixels outside the overlap
        tempzm(ones_map==0) = 0;
        
        normOftemp = norm(tempzm(:));
        % Save current correlation coefficient
        rho = dot(temp(:),wim(:)) / normOftemp / norm(wim(:));
        
        if (i == noi) % the algorithm is executed (noi-1) times
            break;
        end
        
        
        con=cond(C);
        if con>1.0e+15
            disp('->ECC-IC Warning: Badly conditioned Hessian matrix. Check the initialization or the overlap of images.')
        end
        %i_C = inv(C);
        
        % Compute projections of images into G
        Gt = G' * tempzm(:);
        Gw = G' * wim(:);
        
        
        %% ECC closed form solution
        
        vector = C\Gt; % inv(C)*Gt
        
        % Compute lambda parameter
        num = (normOftemp^2 - dot(Gt,vector));
        den = (dot(tempzm(:),wim(:)) - dot(Gw,vector));
        
        lambda = num / den;
        
        % Compute error vector
        imerror = lambda * wim - tempzm;
        
        % Compute the projection of error vector into Jacobian G
        Ge = G' * imerror(:);
        
        % Compute the optimum parameter correction vector
        delta_p = C\Ge; % inv(C)*Ge
        
        if (sum(isnan(delta_p)))>0 %Hessian is close to singular
            disp([' -> ECC-IC algorithms stopped at ' num2str(i) '-th iteration of ' num2str(nol) '-th level due to bad condition of Hessian matrix.']);
            break_flag=1;
            break;
        end
        
        % Update parmaters
        warp = iat_warp_updateIC(warp, delta_p, transform);

    end
    
    if break_flag==1
        break;
    end
    
    % modify the parameteres appropriately for next pyramid level
    if (nol>1)&&(break_flag==0)
        warp = iat_warp_updown(warp, transform,1);
    end
    if levels>1
    fprintf('Done\n');
    end
end

if break_flag==1 % this conditional part is only executed when algorithm stops due to Hessian singularity
    for jj=1:nol-1
        warp = iat_warp_updown(warp, transform,1);
    end
end

%%%