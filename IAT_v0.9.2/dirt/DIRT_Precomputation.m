function pre = DIRT_Precomputation(sI, varargin)
% pre = DIRT_Precomputation(sI, ...)
%
% Do some precomputations for image alignment using the DIRT_Registration
% function.
%
% Mandatory inputs:
%  - sI [image]
%       The source image
%
% Optional inputs (properties):
%  - 'ROI' [struct from the DIRT_Mask2ROI function] Default: user is asked to enter a polygonal ROI  
%       A structure describing the region of interest
%  - 'FM_RC' [m x 2 matrix] Default: []
%       A matrix containing the coordinates of points to be matched to the
%       target image. This makes DIRT enters the FM mode. Only one of ROI and FM_RC can
%       be specified at the same time
%  - 'HalfMaskSize' [scalar] Default: 5
%       * FM mode only *
%       Half the side length of the ROI around each point in FM mode
%  - 'gopt' [struct] Default: []
%       * ROI mode only *
%       A structure parameter to be passed to the DIRT_GModel_* functions
%  - 'no_goptim' [boolean]
%       Toggle optimization of the geometric registration off
%  - 'gmodel' [string] Default: 'Homography' in ROI mode and '3PtHomography' in FM mode
%       The geometric registration model
%       For '3PtHomography', a source and a target point must be given in
%       gopt as column 2-vectors gopt.sq and gopt.tq
%  - 'no_poptim' [boolean]
%       Toggle optimization of the photometric registration off
%  - 'pmodel' [string] Default: 'GainAndBias' for gray-level images
%       The photometric registration model
%  - 'verb' [scalar] Default: 1
%       The verbosity level
%  - 'debug' [boolean]
%       Toggle debug (slow) mode
%  - 'debug2' [boolean]
%       Toggle debug2 (very slow) mode - keeps everything
%  - 'filter_image_der_sigma' [scalar] Default: 1
%       Scale of the Gaussian derivative filter used to differentiate the
%       source image sI
%  - 'no_pnorm' [boolean]
%       Toggle normalization of the pixel values off
%  - 'no_gnorm' [boolean]
%       Toggle normalization of the pixel coordinates off
%
% Outputs:
%  - pre [struct]
%       The sought-after structure, always containing:
%           * optPre [struct]
%               The precomputation options, defined from the properties
%           * sI_roi [image]
%               The source image pixels for the region of interest
%           * sIc_roi and sIr_roi [images]
%               The source image partial derivatives for the region of
%               interest
%           * R_roi and C_roi [columns vectors]
%               The coordinates of the pixels in the region of interest
%               In FM mode, these are cells with size m
%           * SD_roi [matrix]
%               The steepest descent images for the region of interest
%           * invH [matrix]
%               The inverse of the Hessian matrix for the normal equations
%       and that may contain, in debug mode:
%           * sI [image]
%               The source image
%           * sIc and sIr [images]
%               The source image partial derivatives
%           * R_all and C_all [columns vectors]
%               The coordinates of all pixels in the source image
%           * SD_all [matrix]
%               The steepest descent images for all pixels in the source
%               image
%           * H [matrix]
%               The Hessian matrix for the normal equations
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

% check number of parameters
if nargin < 2, error('not enough input arguments'); end;

% setup default parameters and parse input parameters
pre.optPre = ParseVarargin(ABIsRGB(sI), sI, varargin{:});

% print a summary of the input parameters if required
if pre.optPre.verb, DIRT_optPre_PrintSummary(pre.optPre); end;

% image size
nr = size(sI, 1);
nc = size(sI, 2);

% store the source image data
if ~pre.optPre.color
    pre.sI_roi = sI(pre.optPre.ROI.Rind);
else
    sI_1 = sI(:,:,1);
    pre.sI_roi(:, 1) = sI_1(pre.optPre.ROI.Rind);
    sI_2 = sI(:,:,2);
    pre.sI_roi(:, 2) = sI_2(pre.optPre.ROI.Rind);
    sI_3 = sI(:,:,3);
    pre.sI_roi(:, 3) = sI_3(pre.optPre.ROI.Rind);
end;
if pre.optPre.debug
    pre.sI = sI;
end;
if strcmp(pre.optPre.Mode,'FM'), pre.sI = sI; end;

% photometric range normalization (pixel values)
pre.optPre.pnorm_params = DIRT_PNorm_Params(pre.optPre.pnorm, pre.sI_roi);

% photometrically normalized sI_roi
pre.sI_roi = DIRT_PNorm_NormImage(pre.optPre.pnorm_params, pre.sI_roi);

% the image derivatives
if pre.optPre.filter_image_der_sigma == 0
    filter = ABFilters1D('der_x');
else
    filter = ABFilters1D('gauss_der_x', pre.optPre.filter_image_der_sigma);
end;
if ~pre.optPre.color
    sIc = conv2(sI, - filter, 'same');
    sIr = conv2(sI, - filter', 'same');
    pre.sIc_roi = sIc(pre.optPre.ROI.Rind);
    pre.sIr_roi = sIr(pre.optPre.ROI.Rind);
else
    for c = 1:3
        sIc_c = conv2(sI(:, :, c), - filter, 'same');
        pre.sIc_roi(:, c) = sIc_c(pre.optPre.ROI.Rind);
        sIr_c = conv2(sI(:, :, c), - filter', 'same');
        pre.sIr_roi(:, c) = sIr_c(pre.optPre.ROI.Rind);
    end;
end;
if pre.optPre.debug2
    if ~pre.optPre.color
        pre.sIc = sIc;
        pre.sIr = sIr;
    else
        for c = 1:3
            pre.sIc(:, :, c) = conv2(sI(:, :, c), - filter, 'same');
            pre.sIr(:, :, c) = conv2(sI(:, :, c), - filter', 'same');
        end;
    end;
end;

% the image gradient must be normalized, assuming that the normalizing
% transformation can be applied after derivation
pre.sIc_roi = DIRT_PNorm_NormImage(pre.optPre.pnorm_params, pre.sIc_roi);
pre.sIr_roi = DIRT_PNorm_NormImage(pre.optPre.pnorm_params, pre.sIr_roi);

% columns and rows for all pixels in the ROI and in the image
[C_all, R_all] = meshgrid(1:nc, 1:nr);

% columns and rows for pixels of interest
switch pre.optPre.Mode
    case 'ROI'
        pre.C_roi = C_all(pre.optPre.ROI.Rind);
        pre.R_roi = R_all(pre.optPre.ROI.Rind);
        if pre.optPre.debug
            pre.C_all = C_all(:);
            pre.R_all = R_all(:);
        end;
    case 'FM'
        for sj = 1:pre.optPre.m
            if pre.optPre.verb, fprintf('sj = %04d / %04d\n', sj, pre.optPre.m); end;    
            ROI = DIRT_Mask2ROI(DIRT_MaskRectangle(sI, round(pre.optPre.sRC(sj,1)), round(pre.optPre.sRC(sj,2)), pre.optPre.HalfMaskSize, pre.optPre.HalfMaskSize));
            pre.ROI_Rind{sj} = ROI.Rind;
            pre.C_roi{sj} = C_all(ROI.Rind);
            pre.R_roi{sj} = R_all(ROI.Rind);
        end;
    otherwise
        error('unknown mode');
end;

% computation of the Jacobian and Hessian matrices
switch pre.optPre.Mode 
    case 'ROI'
        if pre.optPre.debug2
            [pre.SD_roi, pre.invH, pre.H, pre.SD_all] = DIRT_Precomputation_ROI(pre);
        elseif pre.optPre.debug
            [pre.SD_roi, pre.invH, pre.H] = DIRT_Precomputation_ROI(pre);
        else
            [pre.SD_roi, pre.invH] = DIRT_Precomputation_ROI(pre);
        end;
    case 'FM'
        for j = 1:pre.optPre.m
            pre.optPre.gopt.sq = pre.optPre.sRC(j,:)';
            if pre.optPre.debug
                [pre.SD_roi{j}, pre.invH{j}, pre.H{j}] = DIRT_Precomputation_FM(pre, pre.R_roi{j}, pre.C_roi{j}, pre.ROI_Rind{j});
            else
                [pre.SD_roi{j}, pre.invH{j}] = DIRT_Precomputation_FM(pre, pre.R_roi{j}, pre.C_roi{j}, pre.ROI_Rind{j});
            end;
        end;
    otherwise
        error('unknown mode');
end;

%---
%- Property parsing and default option values
%---
function optPre = ParseVarargin(color, sI, varargin)

% default values
optPre.gopt = [];
optPre.color = color;
optPre.goptim = true;
optPre.poptim = true;
if ~optPre.color
    optPre.pmodel = 'GainAndBias';
else
    optPre.pmodel = 'SingleGainAndBias';
end;
optPre.np = DIRT_PModel_N(optPre.pmodel, optPre.color);
optPre.verb = 1;
optPre.debug = false;
optPre.debug2 = false;
optPre.filter_image_der_sigma = 1;
optPre.gnorm = true;
optPre.pnorm = true;
optPre.HalfMaskSize = 5;
optPre.Mode = 'UnassignedYet';

% varargin parsing
z = 1;
while z <= length(varargin)
    switch(lower(varargin{z}))
        case 'roi'
            if isfield(optPre,'FM_RC'), error('ROI is specified while FM_RC has already been'); end;
            optPre.ROI = varargin{z+1};
            optPre.Mode = 'ROI';
            z = z + 2;
        case 'fm_rc'
            if isfield(optPre,'ROI'), error('FM_RC is specified while ROI has already been'); end;
            optPre.sRC = varargin{z+1};
            optPre.m = size(optPre.sRC, 1);
            optPre.Mode = 'FM';
            optPre.ROI = DIRT_Mask2ROI(DIRT_MaskBorder(0,size(sI,1),size(sI,2)));
            z = z + 2;
        case 'halfmasksize'
            optPre.HalfMaskSize = varargin{z+1};
            z = z + 2;
        case 'gopt'
            optPre.gopt = varargin{z+1};
            z = z + 2;
        case 'no_goptim'
            optPre.goptim = false;
            z = z + 1;
        case 'gmodel'
            optPre.gmodel = varargin{z+1};
            optPre.ng = DIRT_GModel_N(optPre.gmodel);
            z = z + 2;
        case 'no_poptim'
            optPre.poptim = false;
            z = z + 1;
        case 'pmodel'
            optPre.pmodel = varargin{z+1};
            optPre.np = DIRT_PModel_N(optPre.pmodel, optPre.color);
            z = z + 2;
        case 'verb'
            optPre.verb = varargin{z+1};
            z = z + 2;
        case 'debug'
            optPre.debug = true;
            z = z + 1;
        case 'debug2'
            optPre.debug = true;
            optPre.debug2 = true;
            z = z + 1;
        case 'filter_image_der_sigma'
            optPre.filter_image_der_sigma = varargin{z+1};
            z = z + 2;
        case 'no_gnorm'
            optPre.gnorm = false;
            z = z + 1;
        case 'no_pnorm'
            optPre.pnorm = false;
            z = z + 1;
        otherwise
            error(['unknown property ' varargin{z}]);
    end;
end;

% ROI
if strcmp(optPre.Mode,'UnassignedYet') %~isfield(optPre,'ROI') && ~isfield(optPre,'FM_RC')
    fprintf('Enter a polygonal ROI (Region Of Interest)\n');
    optPre.ROI = DIRT_Mask2ROI(DIRT_MaskManual(sI)); 
    optPre.Mode = 'ROI';
end;
% gmodel
if ~isfield(optPre,'gmodel')
    switch optPre.Mode
        case 'ROI'
            optPre.gmodel = 'Homography';
        case 'FM'
            optPre.gmodel = '3PtHomography';
    end;
else
    if ~strcmp(optPre.gmodel,'Homography') ...
        && ~strcmp(optPre.gmodel,'3PtHomography') ...
        && ~strcmp(optPre.gmodel,'Affine') ...
        && ~strcmp(optPre.gmodel,'Rt')
        error('Unknown gmodel');
    end;
end;
optPre.ng = DIRT_GModel_N(optPre.gmodel);
% finish initialization
% number of parameters to be optimized
if optPre.goptim, optPre.no = optPre.ng; else optPre.no = 0; end;
if optPre.poptim, optPre.no = optPre.no + optPre.np; end;

%---
%- Computation of the constant Jacobian and Hessian for some pixels of interest in ROI mode
%---
function [SD_roi, invH, oH, oSD_all] = DIRT_Precomputation_ROI(pre)

% the steepest descent images for the geometric registration parameters
if pre.optPre.goptim    
    [Gr, Gc] = DIRT_GModel_GradientVectors(pre.optPre.gmodel, pre.R_roi, pre.C_roi, pre.optPre.gopt);

    if ~pre.optPre.color    
        % in the grey-level case, SD_roi is (Rn x ng)
        SD_roi = repmat(pre.sIr_roi, 1, pre.optPre.ng).*Gr + repmat(pre.sIc_roi, 1, pre.optPre.ng).*Gc;
    else
        % in the color case, the SD_roi interlaces the RGB channels and is
        % (3*Rn x ng)
        SD_roi = zeros(pre.optPre.ROI.Rn*3, pre.optPre.ng);
        for c = 1:3
            SD_roi(c:3:end, :) = repmat(pre.sIr_roi(:, c), 1, pre.optPre.ng).*Gr + repmat(pre.sIc_roi(:, c), 1, pre.optPre.ng).*Gc;
        end;
    end;
    
    if nargout > 3
        [Gr, Gc] = DIRT_GModel_GradientVectors(pre.optPre.gmodel, pre.R_all, pre.C_all, pre.optPre.gopt);

        if ~pre.optPre.color
            % in the grey-level case, the columns of SD_all are the
            % vectorized steepest descent images
            SD_all = repmat(pre.sIr(:), 1, pre.optPre.ng).*Gr + repmat(pre.sIc(:), 1, pre.optPre.ng).*Gc;
        else
            % in the color case, a third dimension is added to the matrix,
            % so that the steepest descent images for all three channels
            % are represented. SD_all is thus (3*Rn x ng x 3)
            sIr_1 = pre.sIr(:,:,1);
            sIc_1 = pre.sIc(:,:,1);
            SD_all(:,:,1) = repmat(sIr_1(:), 1, pre.optPre.ng).*Gr + repmat(sIc_1(:), 1, pre.optPre.ng).*Gc;
            sIr_2 = pre.sIr(:,:,2);
            sIc_2 = pre.sIc(:,:,2);
            SD_all(:,:,2) = repmat(sIr_2(:), 1, pre.optPre.ng).*Gr + repmat(sIc_2(:), 1, pre.optPre.ng).*Gc;
            sIr_3 = pre.sIr(:,:,3);
            sIc_3 = pre.sIc(:,:,3);
            SD_all(:,:,3) = repmat(sIr_3(:), 1, pre.optPre.ng).*Gr + repmat(sIc_3(:), 1, pre.optPre.ng).*Gc;
        end;
    end;
else
    SD_all = [];
end;

% the steepest descent images for the photometric part
if pre.optPre.poptim
    SD_roi = [ SD_roi DIRT_PModel_GradientVectors(pre.optPre.pmodel, pre.optPre.color, pre.sI_roi) ];
    if pre.optPre.debug2
        if ~pre.optPre.color
            SD_all = [ SD_all DIRT_PModel_GradientVectors(pre.optPre.pmodel, pre.optPre.color, pre.sI(:)) ];
        else
            % implementation to be done
%             sI_1 = pre.sI(:,:,1);
%             sI_2 = pre.sI(:,:,2);
%             sI_3 = pre.sI(:,:,3);
%             v_sI(:, 1) = sI_1(:);
%             v_sI(:, 2) = sI_2(:);
%             v_sI(:, 3) = sI_3(:);
        end;
    end;
end;

% the (Gauss-Newton approximation to the) Hessian matrix
if pre.optPre.goptim || pre.optPre.poptim
    H = SD_roi' * SD_roi;
    if cond(H) > 10^15
        error('The source ROI does not seem to contain enough information to register -- try a simpler gmodel/pmodel');
    end;
    invH = inv(H);
end;

if nargout > 2, oH = H; end;
if nargout > 3, oSD_all = SD_all; end;

%---
%- Computation of the constant Jacobian and Hessian for some pixels of interest in FM mode
%---
function [SD_roi, invH, oH] = DIRT_Precomputation_FM(pre, R_roi, C_roi, Rind)

Rn = length(Rind);

% the steepest descent images for the geometric registration parameters
if pre.optPre.goptim    
    [Gr, Gc] = DIRT_GModel_GradientVectors(pre.optPre.gmodel, R_roi, C_roi, pre.optPre.gopt);

    if ~pre.optPre.color    
        % in the grey-level case, SD_roi is (Rn x ng)
        SD_roi = repmat(pre.sIr_roi(Rind), 1, pre.optPre.ng).*Gr + repmat(pre.sIc_roi(Rind), 1, pre.optPre.ng).*Gc;
    else
        % in the color case, the SD_roi interlaces the RGB channels and is
        % (3*Rn x ng)
        SD_roi = zeros(Rn*3, pre.optPre.ng);
        for c = 1:3
            SD_roi(c:3:end, :) = repmat(pre.sIr_roi(Rind,c), 1, pre.optPre.ng).*Gr + repmat(pre.sIc_roi(Rind,c), 1, pre.optPre.ng).*Gc;
        end;
    end;
end;
    
% the steepest descent images for the photometric part
if pre.optPre.poptim
    SD_roi = [ SD_roi DIRT_PModel_GradientVectors(pre.optPre.pmodel, pre.optPre.color, pre.sI_roi(Rind,:)) ];
end;

% the (Gauss-Newton approximation to the) Hessian matrix
if pre.optPre.goptim || pre.optPre.poptim
    H = SD_roi' * SD_roi;
    invH = inv(H);
end;

if nargout > 2, oH = H; end;
