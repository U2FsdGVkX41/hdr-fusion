function reg = DIRT_Registration(tI, pre, varargin)
% reg = DIRT_Registration(tI, pre, ...)
%
% Register the target image tI to the source image sI in pre, using the 
% pixels in some region of interest in pre. The registration model is:
%   sI[q] = P(tI[G(q)]) + noise
% where P is a photometric transformation and G a geometric one.
% The function iteratively minimizes the least squares discrepancy in the
% above model.
%
% Mandatory inputs:
%  - tI [image]
%       The target image tI to be registered
%  - pre [struct from DIRT_Precomputation]  
%       A structure containing precomputed quantities, obtained from the
%       DIRT_Precomputation function (it contains the source image sI and 
%       the region of interest)
%
% Optional inputs (properties):
%  - 'FM_RC' [m x 2 matrix] Default: []
%       * FM mode only *
%       A matrix containing the coordinates of points to be matched.
%  - 'gopt' Default: []
%       * ROI mode only *
%       A structure parameter to be passed to the DIRT_GModel_* functions
%       Note that it is added to the pre.gopt structure
%  - 'ginit' [vector or matrix] Default: identity
%       The initial geometric registration
%  - 'pinit' [vector or matrix] Default: identity
%       The initial photometric registration
%  - 'verb' [scalar] Default: 1 in ROI mode and 0 in FM mode
%       The verbosity level
%  - 'FM_verb' [scalar] Default: 1
%       * FM mode only *
%       The verbosity level for FM messages
%  - 'min_step_size' [scalar] Default: 1e-8
%       Minimum size of the update vector for iterating
%  - 'max_nb_it' [scalar] Default: 100
%       Maximum number of iterations
%  - 'pause_after_update' [boolean] Default: false
%       Pause after each parameter update
%  - 'unfilled_pixels_diff' [scalar or vector] Default: 0 or [ 0 0 0 ]
%       Default color for the pixels which are not filled in the difference
%       image. Only makes sense if pre.optPre.debug = true
%  - 'unfilled_pixels_warped' [scalar or vector] Default: 0 or [ 0 0 0 ]
%       Default color for the pixels which are not filled in the warped
%       image. Only makes sense if pre.optPre.debug = true
%  - 'write_res_path' [string] Default: []
%       Write the results to the specified path
%  - 'no_max_nb_it' [boolean]
%       Do not stop after max_nb_it is passed
%  - 'no_min_step_size' [boolean]
%       Do not stop after min_step_size is passed
%  - 'max_init_error' [scalar] Default: +inf
%       Stop if the initial error is larger than max_init_error
%       Useful mainly in FM mode
%  - 'max_disparity' [scalar] Default: +inf
%       * FM mode only *
%       Do not compute OptSSD if the disparity between two points is larger
%       than this threshold
%       
% Outputs:
%  - reg [struct]
%   In ROI mode:
%       Structure always containing:
%           * err [boolean]
%               Error in the registration process
%           * g [vector or matrix]
%               Geometric registration parameters
%           * p [vector or matrix]
%               Photometric registration parameters
%           * it [scalar]
%               Number of iterations
%           * cpu_time [scalar]
%               Total cpu time in seconds
%           * e_roi [scalar]
%               RMS fitting error for the ROI (in pixel value units)
%           * e_roi_init [scalar]
%               Initial RMS fitting error for the ROI (in pixel value units)
%       Structure that may contain:
%           * allReg [list of struct]
%               This parameter is returned if 'debug' has been activated in
%               DIRT_Precomputation
%               allReg{t} with t = 1 ... reg.it contains information about 
%               the registration through at iteration t, in debug mode:
%                   - g, p, cpu_time, e_roi [see above]
%                   - e [scalar]
%                       RMS fitting error for the whole image (in pixel
%                       value units)
%               and in debug2 mode:
%                   - wI [image]
%                       Warped image wI[q] = P(tI[G(q)])
%                   - dI [image]
%                       Difference image dI = wI - sI
%   In FM mode:
%       * DispMat [matrix]
%           The disparity matrix (DispMat(sj,tj) is the distance between
%           the sj-th and tj-th source and target points)
%       * RawSSD [matrix]
%           The raw SSD matrix
%       * OptSSD [matrix]
%           The optimized SSD matrix
%       * nb_OptSSD_comp [scalar]
%           Number of optimized SSD computed (i.e. number of point pairs
%           which passed the maximum disparity and maximum initial error
%           tests)
%                       
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

% check that the number of channels in the source and target images are identical
if ABIsRGB(tI) ~= pre.optPre.color, error('the number of channels of the source and target images is different'); end;

% normalization
tI_norm = DIRT_PNorm_NormImage(pre.optPre.pnorm_params, tI);

% check number of parameters
if nargin < 2, error('not enough input arguments'); end;

% setup default parameters and parse input parameters
optReg = ParseVarargin(pre.optPre, varargin{:});

% run registration
switch pre.optPre.Mode
    case 'ROI'
        reg = DIRT_Registration_ROI(pre, optReg, tI_norm);
    case 'FM'
        reg.nb_OptSSD_comp = 0;
        
        % make the disparity matrix
        reg.DispMat = sqrt( ...
            (repmat(pre.optPre.sRC(:,1),1,optReg.m) - repmat(optReg.tRC(:,1)',pre.optPre.m,1)).^2 + ...
            (repmat(pre.optPre.sRC(:,2),1,optReg.m) - repmat(optReg.tRC(:,2)',pre.optPre.m,1)).^2 );
        for sj = 1:pre.optPre.m
            if optReg.FM_verb, fprintf('sj = %04d / %04d\n', sj, pre.optPre.m); end;    
            optReg.gopt.sq = pre.optPre.sRC(sj,:)';
            for tj = 1:optReg.m
%                if optReg.FM_verb, fprintf('  tj = %04d / %04d\n', tj, optReg.m); end;

                r1 = round(pre.optPre.sRC(sj,1));
                c1 = round(pre.optPre.sRC(sj,2));
                r2 = round(optReg.tRC(tj,1));
                c2 = round(optReg.tRC(tj,2));
                hs = pre.optPre.HalfMaskSize;
                
                % compute pixel level SSD
                if ~pre.optPre.color
                    SSD = norm(pre.sI(r1-hs:r1+hs,c1-hs:c1+hs) - tI_norm(r2-hs:r2+hs,c2-hs:c2+hs),'fro') / (2*hs+1);
                else
                    SSD = sqrt( ...
                        norm(pre.sI(r1-hs:r1+hs,c1-hs:c1+hs,1) - tI_norm(r2-hs:r2+hs,c2-hs:c2+hs,1),'fro')^2 + ...
                        norm(pre.sI(r1-hs:r1+hs,c1-hs:c1+hs,2) - tI_norm(r2-hs:r2+hs,c2-hs:c2+hs,2),'fro')^2 + ...
                        norm(pre.sI(r1-hs:r1+hs,c1-hs:c1+hs,3) - tI_norm(r2-hs:r2+hs,c2-hs:c2+hs,3),'fro') ) / (3*(2*hs+1));
                end;

                % disparity and pixel level SSD check
                if reg.DispMat(sj,tj) < optReg.max_disparity && SSD < optReg.max_init_error
                    reg.nb_OptSSD_comp = reg.nb_OptSSD_comp + 1;
                    optReg.gopt.tq = optReg.tRC(tj,:)';
                    optReg.ginit = DIRT_GModel_Identity(pre.optPre.gmodel, optReg.gopt);
                    tmp_reg = DIRT_Registration_FM(pre, optReg, tI_norm, sj);
                    if ~tmp_reg.err
                        reg.RawSSD(sj,tj) = tmp_reg.e_roi_init;
                        reg.OptSSD(sj,tj) = tmp_reg.e_roi;
                    else
                        reg.RawSSD(sj,tj) = SSD;
                        reg.OptSSD(sj,tj) = SSD;
                    end;
                else
                    reg.RawSSD(sj,tj) = SSD;
                    reg.OptSSD(sj,tj) = SSD;
                end;
            end;
        end;
    otherwise
        error('unknown mode');
end;

%---
%- Property parsing and default option values
%---
function optReg = ParseVarargin(optPre, varargin)

% default values
optReg.gopt = optPre.gopt;
optReg.FM_verb = 1;
optReg.min_step_size = 1e-8;
optReg.pinit = DIRT_PModel_Identity(optPre.pmodel, optPre.color);
optReg.max_nb_it = 100;
optReg.pause_after_update = false;
if ~optPre.color
    optReg.unfilled_pixels_warped = 0;
    optReg.unfilled_pixels_diff = 0;
else
    optReg.unfilled_pixels_warped = [0 0 0];
    optReg.unfilled_pixels_diff = [0 0 0];
end;
optReg.write_res_path = [];
optReg.no_max_nb_it = false;
optReg.no_min_step_size = false;
optReg.max_init_error = inf;
optReg.max_disparity = inf;

% varargin parsing
z = 1;
while z <= length(varargin)
    switch(lower(varargin{z}))
        case 'fm_rc'
            if strcmp(optPre.Mode,'ROI'), error('the rc property is useless in ROI mode'); end;
            optReg.tRC = varargin{z+1};
            optReg.m = size(optReg.tRC,1);
            z = z + 2;
        case 'gopt'
            gopt = varargin{z+1};
            z = z + 2;
            for f = fieldnames(gopt)
                tmp = strcat('optReg.gopt.',f,' = gopt.',f,';');
                eval(tmp{1});
            end;
        case 'verb'
            optReg.verb = varargin{z+1};
            z = z + 2;
        case 'FM_verb'
            optReg.FM_verb = varargin{z+1};
            z = z + 2;
        case 'min_step_size'
            optReg.min_step_size = varargin{z+1};
            z = z + 2;
        case 'ginit'
            optReg.ginit = varargin{z+1};
            z = z + 2;
        case 'pinit'
            optReg.pinit = varargin{z+1};
            z = z + 2;
        case 'max_nb_it'
            optReg.max_nb_it = varargin{z+1};
            z = z + 2;
        case 'pause_after_update'
            optReg.pause_after_update = true;
            z = z + 1;
        case 'unfilled_pixels_warped'
            optReg.unfilled_pixels_warped = varargin{z+1};
            z = z + 2;
        case 'unfilled_pixels_diff'
            optReg.unfilled_pixels_diff = varargin{z+1};
            z = z + 2;            
        case 'write_res_path'
            optReg.write_res_path = varargin{z+1};
            z = z + 2;            
        case 'no_max_nb_it'
            optReg.no_max_nb_it = true;
            z = z + 1;
        case 'no_min_step_size'
            optReg.no_min_step_size = true;
            z = z + 1;
        case 'max_init_error'
            optReg.max_init_error = varargin{z+1};
            z = z + 2;
        case 'max_disparity'
            optReg.max_disparity = varargin{z+1};
            z = z + 2;
        otherwise
            error(['unknown property ' varargin{z}]);
    end;
end;

if ~isfield(optReg, 'ginit')
    if strcmp(optPre.Mode,'ROI')
        optReg.ginit = DIRT_GModel_Identity(optPre.gmodel, optReg.gopt);
    end;
end;

if strcmp(optPre.Mode,'FM') && ~isfield(optReg,'tRC')
    error('the FM_RC property is required in FM mode');
end;

if ~isfield(optReg,'verb')
    optReg.verb = strcmp(optPre.Mode,'ROI');
end;

%---
%- Main registration subfunction in ROI mode
%---
function reg = DIRT_Registration_ROI(pre, optReg, tI_norm)

% make directories for writing images
if ~isempty(optReg.write_res_path)
    % base path
    if pre.optPre.poptim, pmodelname = pre.optPre.pmodel; else pmodelname = [ pre.optPre.pmodel 'NoOpt' ]; end;
    if pre.optPre.goptim, gmodelname = pre.optPre.gmodel; else gmodelname = [ pre.optPre.gmodel 'NoOpt' ]; end;        
    mkdir(optReg.write_res_path, [ gmodelname '_' pmodelname sprintf('_%08d', length(pre.sI_roi)) ]);
    optReg.write_basepath = [ optReg.write_res_path '/' gmodelname '_' pmodelname sprintf('_%08d', length(pre.sI_roi)) ];
    if pre.optPre.debug
        % difference images
        mkdir(optReg.write_basepath, 'dI');
        optReg.write_dIfmtfn = [ optReg.write_basepath '/dI/%03d.png' ];
        % warped images
        mkdir(optReg.write_basepath, 'wI');
        optReg.write_wIfmtfn = [ optReg.write_basepath '/wI/%03d.png' ];
        % warped images without photometric model
        mkdir(optReg.write_basepath, 'wI_no_pmodel');
        optReg.write_wInopmodelfmtfn = [ optReg.write_basepath '/wI_no_pmodel/%03d.png' ];
        % warped images without geometric model
        mkdir(optReg.write_basepath, 'wI_no_gmodel');
        optReg.write_wInogmodelfmtfn = [ optReg.write_basepath '/wI_no_gmodel/%03d.png' ];
    end;        
end;

% initialize the registration
reg.g = optReg.ginit;
reg.p = DIRT_PNorm_Norm_PModel(pre.optPre.pnorm_params, pre.optPre.pmodel, pre.optPre.color, optReg.pinit);

% iteration counter
reg.it = 0;

% total computational time
reg.cpu_time = 0;

% setup text display
if optReg.verb || pre.optPre.debug
    fprintf('Iteration | RMS error (ROI) | cpu time (s)')
end;
if pre.optPre.debug
    fprintf(' | RMS error (total)');
end;
if optReg.verb || pre.optPre.debug
    DIRT_GModel_DispInit(pre.optPre.gmodel);
    DIRT_PModel_DispInit(pre.optPre.pmodel, pre.optPre.color);
    fprintf('\n');
end;

% structure for storing the registration through the iterations
if pre.optPre.debug
    reg.allReg = {};
end;

% main loop
converged = false;
while true    
    
    % compute the warped and difference images and the RMS error
    tic;
    [reg.err, wI_roi, dI_roi, reg.e_roi] = DIRT_Warp(pre.optPre.gmodel, reg.g, pre.optPre.pmodel, reg.p, ...
        pre.optPre.color, pre.R_roi, pre.C_roi, tI_norm, false, 0, pre.sI_roi);
    reg.cpu_time = reg.cpu_time + toc;
    if reg.err
        fprintf('DIRT:Warping attempted outside the target image, exiting\n');
        fprintf('DIRT:You may try to improve initialization with DIRT_GModel_Init and/or shrink the ROI with DIRT_MaskBorder\n');
        return; 
    end;

    % de-normalize the RMS error
    reg.e_roi = DIRT_PNorm_UnNormRMSError(pre.optPre.pnorm_params, reg.e_roi);

    % store initial error and check it
    if ~reg.it
        reg.e_roi_init = reg.e_roi;
        if reg.e_roi_init >= optReg.max_init_error
            if optReg.verb, fprintf('Initial error too large\n'); end;
            break;
        end;
    end;

    % de-normalize the photometric registration parameters
    p_unnorm = DIRT_PNorm_UnNorm_PModel(pre.optPre.pnorm_params, pre.optPre.pmodel, pre.optPre.color, reg.p);
    
    if pre.optPre.debug
        [err, wI, dI, e, wI_no_pmodel, wI_no_gmodel] = DIRT_Warp(pre.optPre.gmodel, reg.g, ...,
            pre.optPre.pmodel, p_unnorm, pre.optPre.color, pre.R_all, pre.C_all, tI_norm, true, ...,
            optReg.unfilled_pixels_warped, pre.sI, optReg.unfilled_pixels_diff);
    end;

    % store the current estimate
    if pre.optPre.debug
        reg.allReg{end+1}.g = reg.g;
        reg.allReg{end}.p = p_unnorm;
        reg.allReg{end}.e = e;
        reg.allReg{end}.cpu_time = reg.cpu_time;
        reg.allReg{end}.e_roi = reg.e_roi;
        if pre.optPre.debug2
            reg.allReg{end}.wI = wI;
            reg.allReg{end}.dI = dI;
        end;        
    end;
    
    % write images to disk
    if ~isempty(optReg.write_res_path) && pre.optPre.debug
        imwrite(uint8(abs(dI)), sprintf(optReg.write_dIfmtfn, reg.it)); 
        imwrite(uint8(abs(wI)), sprintf(optReg.write_wIfmtfn, reg.it)); 
        imwrite(uint8(abs(wI_no_pmodel)), sprintf(optReg.write_wInopmodelfmtfn, reg.it)); 
        imwrite(uint8(abs(wI_no_gmodel)), sprintf(optReg.write_wInogmodelfmtfn, reg.it)); 
    end;
    
    % display information
    if optReg.verb || pre.optPre.debug
        fprintf('%03d       | %012.7f    | %012.7f', reg.it, reg.e_roi, reg.cpu_time);
    end;
    if pre.optPre.debug
        fprintf(' | %012.7f     ', e);
    end;
    if optReg.verb || pre.optPre.debug
        DIRT_GModel_Disp(pre.optPre.gmodel, reg.g);
        DIRT_PModel_Disp(pre.optPre.pmodel, pre.optPre.color, p_unnorm);
        fprintf('\n');
    end;
    
    % if no optimization required
    if ~pre.optPre.goptim && ~pre.optPre.poptim 
        if optReg.verb, fprintf('No optimization required\n'); end;
        converged = true; 
    end;
    
    % convergence test
    if converged, break; end;
    
    tic;

    % compute the right hand side of the normal equations
    b = sum(pre.SD_roi' .* repmat(dI_roi', pre.optPre.no, 1), 2);

    % compute the update vector
    d = pre.invH * b;

    % norm of the update vector for convergence test
    normd = norm(d);

    % update the current geometric estimate
    if pre.optPre.goptim
        dg = d(1:pre.optPre.ng);
        d = d(pre.optPre.ng+1:end);
        reg.g = DIRT_GModel_Update(pre.optPre.gmodel, reg.g, dg, optReg.gopt);
    end;

    % update the current photometric estimate
    if pre.optPre.poptim
        dp = d(1:pre.optPre.np);
        reg.p = DIRT_PModel_Update(pre.optPre.pmodel, pre.optPre.color, reg.p, dp);
    end;
    
    reg.cpu_time = reg.cpu_time + toc;    
    
    % convergence test
    if normd < optReg.min_step_size 
        if optReg.verb, fprintf('Minimum step size reached\n'); end;            
        converged = ~optReg.no_min_step_size; 
    end;

    % increase iteration counter
    reg.it = reg.it + 1;
    
    % maximum number of iterations reached
    if reg.it >= optReg.max_nb_it 
        if optReg.verb, fprintf('Maximum number of iterations reached\n'); end;
        converged = ~optReg.no_max_nb_it;
    end;
    
    % optional pause for debugging purposes
    if optReg.pause_after_update, pause; end;

end; % main loop

% de-normalize the photometric registration parameters
reg.p = DIRT_PNorm_UnNorm_PModel(pre.optPre.pnorm_params, pre.optPre.pmodel, pre.optPre.color, reg.p);

% save the registration structure and options, and pre-computations
if ~isempty(optReg.write_res_path)
    save([optReg.write_basepath '/reg_optReg_pre.mat'], 'reg', 'optReg', 'pre');
end;

%---
%- Main registration subfunction in FM mode
%---
function reg = DIRT_Registration_FM(pre, optReg, tI_norm, sj)

% % make directories for writing images
% if ~isempty(optReg.write_res_path)
%     % base path
%     if pre.optPre.poptim, pmodelname = pre.optPre.pmodel; else pmodelname = [ pre.optPre.pmodel 'NoOpt' ]; end;
%     if pre.optPre.goptim, gmodelname = pre.optPre.gmodel; else gmodelname = [ pre.optPre.gmodel 'NoOpt' ]; end;        
%     mkdir(optReg.write_res_path, [ gmodelname '_' pmodelname sprintf('_%08d', pre.ROI.Rn) ]);
%     optReg.write_basepath = [ optReg.write_res_path '/' gmodelname '_' pmodelname sprintf('_%08d', pre.ROI.Rn) ];
%     if pre.optPre.debug
%         % difference images
%         mkdir(optReg.write_basepath, 'dI');
%         optReg.write_dIfmtfn = [ optReg.write_basepath '/dI/%03d.png' ];
%         % warped images
%         mkdir(optReg.write_basepath, 'wI');
%         optReg.write_wIfmtfn = [ optReg.write_basepath '/wI/%03d.png' ];
%         % warped images without photometric model
%         mkdir(optReg.write_basepath, 'wI_no_pmodel');
%         optReg.write_wInopmodelfmtfn = [ optReg.write_basepath '/wI_no_pmodel/%03d.png' ];
%         % warped images without geometric model
%         mkdir(optReg.write_basepath, 'wI_no_gmodel');
%         optReg.write_wInogmodelfmtfn = [ optReg.write_basepath '/wI_no_gmodel/%03d.png' ];
%     end;        
% end;

% initialize the registration
reg.g = optReg.ginit;
reg.p = DIRT_PNorm_Norm_PModel(pre.optPre.pnorm_params, pre.optPre.pmodel, pre.optPre.color, optReg.pinit);

% iteration counter
reg.it = 0;

% total computational time
reg.cpu_time = 0;

% setup text display
if optReg.verb || pre.optPre.debug
    fprintf('Iteration | RMS error (ROI) | cpu time (s)')
end;
if pre.optPre.debug
    fprintf(' | RMS error (total)');
end;
if optReg.verb || pre.optPre.debug
    DIRT_GModel_DispInit(pre.optPre.gmodel);
    DIRT_PModel_DispInit(pre.optPre.pmodel, pre.optPre.color);
    fprintf('\n');
end;

% structure for the storing the registration through the iterations
if pre.optPre.debug
    reg.allReg = {};
end;

% main loop
converged = false;
while true    
    
    % compute the warped and difference images and the RMS error
    tic;
    [reg.err, wI_roi, dI_roi, reg.e_roi] = DIRT_Warp(pre.optPre.gmodel, reg.g, pre.optPre.pmodel, reg.p, ...
        pre.optPre.color, pre.R_roi{sj}, pre.C_roi{sj}, tI_norm, false, 0, pre.sI_roi(pre.ROI_Rind{sj},:));
    reg.cpu_time = reg.cpu_time + toc;
    if reg.err
        fprintf('Warping attempted outside the target image, exiting\n');
        fprintf('You may try to improve initialization with DIRT_GModel_Init and/or shrink the ROI with DIRT_MaskBorder\n');
        return; 
    end;
    
    % de-normalize the RMS error
    reg.e_roi = DIRT_PNorm_UnNormRMSError(pre.optPre.pnorm_params, reg.e_roi);

    % store initial error and check it
    if reg.it == 0
        reg.e_roi_init = reg.e_roi;
        if reg.e_roi_init >= optReg.max_init_error
            if optReg.verb, fprintf('Initial error too large\n'); end;
            break;
        end;
    end;
    
    % de-normalize the photometric registration parameters
    p_unnorm = DIRT_PNorm_UnNorm_PModel(pre.optPre.pnorm_params, pre.optPre.pmodel, pre.optPre.color, reg.p);
    
%     if pre.optPre.debug
%         [err, wI, dI, e, wI_no_pmodel, wI_no_gmodel] = DIRT_Warp(pre.optPre.gmodel, reg.g, ...,
%             pre.optPre.pmodel, p_unnorm, pre.optPre.color, pre.R_all, pre.C_all, tI, true, ...,
%             optReg.unfilled_pixels_warped, pre.sI, optReg.unfilled_pixels_diff);
%     end;

%     % store the current estimate
%     if pre.optPre.debug
%         reg.allReg{end+1}.g = reg.g;
%         reg.allReg{end}.p = p_unnorm;
%         reg.allReg{end}.e = e;
%         reg.allReg{end}.cpu_time = reg.cpu_time;
%         reg.allReg{end}.e_roi = reg.e_roi;
%         if pre.optPre.debug2
%             reg.allReg{end}.wI = wI;
%             reg.allReg{end}.dI = dI;
%         end;        
%     end;
    
%     % write images to disk
%     if ~isempty(optReg.write_res_path) && pre.optPre.debug
%         imwrite(uint8(abs(dI)), sprintf(optReg.write_dIfmtfn, reg.it)); 
%         imwrite(uint8(abs(wI)), sprintf(optReg.write_wIfmtfn, reg.it)); 
%         imwrite(uint8(abs(wI_no_pmodel)), sprintf(optReg.write_wInopmodelfmtfn, reg.it)); 
%         imwrite(uint8(abs(wI_no_gmodel)), sprintf(optReg.write_wInogmodelfmtfn, reg.it)); 
%     end;
    
    % display information
    if optReg.verb || pre.optPre.debug
        fprintf('%03d       | %012.7f    | %012.7f', reg.it, reg.e_roi, reg.cpu_time);
    end;
    if pre.optPre.debug
        fprintf(' | %012.7f     ', e);
    end;
    if optReg.verb || pre.optPre.debug
        DIRT_GModel_Disp(pre.optPre.gmodel, reg.g);
        DIRT_PModel_Disp(pre.optPre.pmodel, pre.optPre.color, p_unnorm);
        fprintf('\n');
    end;
    
    % if no optimization required
    if ~pre.optPre.goptim && ~pre.optPre.poptim 
        if optReg.verb, fprintf('No optimization required\n'); end;
        converged = true; 
    end;
    
    % convergence test
    if converged, break; end;
    
    tic;

    % compute the right hand side of the normal equations
    b = sum(pre.SD_roi{sj}' .* repmat(dI_roi', pre.optPre.no, 1), 2);

    % compute the update vector
    d = pre.invH{sj} * b;

    % norm of the update vector for convergence test
    normd = norm(d);

    % update the current geometric estimate
    if pre.optPre.goptim
        dg = d(1:pre.optPre.ng);
        d = d(pre.optPre.ng+1:end);
        reg.g = DIRT_GModel_Update(pre.optPre.gmodel, reg.g, dg, optReg.gopt);
    end;

    % update the current photometric estimate
    if pre.optPre.poptim
        dp = d(1:pre.optPre.np);
        reg.p = DIRT_PModel_Update(pre.optPre.pmodel, pre.optPre.color, reg.p, dp);
    end;
    
    reg.cpu_time = reg.cpu_time + toc;    
    
    % convergence test
    if normd < optReg.min_step_size 
        if optReg.verb, fprintf('Minimum step size reached\n'); end;            
        converged = ~optReg.no_min_step_size; 
    end;

    % increase iteration counter
    reg.it = reg.it + 1;
    
    % maximum number of iterations reached
    if reg.it >= optReg.max_nb_it 
        if optReg.verb, fprintf('Maximum number of iterations reached\n'); end;
        converged = ~optReg.no_max_nb_it;
    end;
    
    % optional pause for debugging purposes
    if optReg.pause_after_update, pause; end;

end; % main loop

% de-normalize the photometric registration parameters
reg.p = DIRT_PNorm_UnNorm_PModel(pre.optPre.pnorm_params, pre.optPre.pmodel, pre.optPre.color, reg.p);

% % save the registration structure and options, and pre-computations
% if ~isempty(optReg.write_res_path)
%     save([optReg.write_basepath '/reg_optReg_pre.mat'], 'reg', 'optReg', 'pre');
% end;
