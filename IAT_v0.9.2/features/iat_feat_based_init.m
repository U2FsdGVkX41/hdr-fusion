function [warp] = iat_feat_based_init(transform, method, varargin)
% iat_feat_based_init.m implements a feature based warp initialization
% according to a specified method.
%
% /Standard input
% transform:            The adopted transform
% method:               The initialization method. Valid values are:
%                       'LS','desc-LS','ext-LS','RANSAC','desc-RANSAC',
%                       'ext-RANSAC'.
% 
% /Variable input
%
% >>for method='LS':
% image:                The initial image
% template:             The template image
% distRatio (optional): Ratio between the first and second best
%                       descriptor match. If their ratio is bigger than 
%                       this value, then the match is not confirmed.
%
% Description:
% The initial warp is determined by extracting descriptors with the
% built-in (SURF) detector, match the descriptors, and obtain the warp with
% the Least-Squares algorithm.
%
%
%
% >>for method='desc-LS':
% iDesc:                Image Descriptors (image_points x descriptors_size)
% tDesc:                Template Descriptors (template_points x descriptor_size)
% iPts:                 Image points (2xN or 3xN, corresponding to iDesc)
% tPts:                 Template points (2xN or 3xN, corresponding to tDesc)
%
% distRatio (optional): Ratio between the first and second best
%                       descriptor match. If their ratio is bigger than 
%                       this value, then the match is not confirmed.
%
% Description:
% The initial warp is determined by external descriptors and points that
% correspond. After the matching of the descriptors, the correspondencies
% are given to the Least-Squares algorithm as input in order to obtain a 
% final estimation.
%
%
%
% >>for method='ext-LS':
% iPts:                 Image points (2xN or 3xN)
% tPts:                 Template points (2xN or 3xN)
%
% Description:
% The initial warp is determined by external correspondencies as input to
% the Least-Squares algorithm.
%
% >>for method='RANSAC':
% image:                The initial image
% template:             The template image
% {RANSAC parameters} (optional)
% distRatio (optional): Ratio between the first and second best descriptor 
%                       match. If their ratio is bigger than this value, 
%                       then the match is not confirmed.
%
% Description:
% The initial warp is determined by extracting descriptors with the
% built-in (SURF) detector, match the descriptors and obtain the warp with
% the RANSAC algorithm.
%
% >>for method='desc-RANSAC':
% iDesc:                Image Descriptors (image_points x descriptors_size)
% tDesc:                Template Descriptors (template_points x descriptor_size)
% iPts:                 Image points (2xN or 3xN, corresponding to iDesc)
% tPts:                 Template points (2xN or 3xN, corresponding to tDesc)
% RANSAC parameters} (optional)
% distRatio (optional): Ratio between the first and second best descriptor 
%                       match. If their ratio is bigger than this value, 
%                       then the match is not confirmed.
%
% Description:
% The initial warp is determined by external descriptors and points that
% correspond. After the matching of the descriptors, the correspondencies
% are given to the RANSAC algorithm as input in order to obtain a final
% estimation.
%
%
%
% >>for method='ext-RANSAC':
% iPts:                 Image points (2xN or 3xN)
% tPts:                 Template points (2xN or 3xN)
% {RANSAC parameters} (optional)
%
% Description:
% The initial warp is determined by external correspondencies as input to
% the RANSAC algorithm.
%
% {RANSAC parameters}:
% rMaxIter:             The maximum number of iterations for RANSAC
% rMaxInvCount:         The maximum number of serial degenerated MSS picks. 
%                       If exceeded, the algorithm is terminated.
% rOutFreeProb:      The outlier-free probability for a consensus set
% rErrorThr:            The error threshold for a data-point to be 
%                       considered as an inlier.
%
% NOTE: Optional parameters are passed as NAME<->VALUE pairs after the
% required parameters, for example:
% H=iat_feat_based_init('affine','RANSAC',img,tmp,'rMaxIter',100,'distRatio',.7);
%
% Output
% -->warp: The extracted warp
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
%% Essential input error checking
if ~(ischar(transform))
    error('iat_feat_based_init.m: Input "transform" must be a string');
end

if ~(ischar(method))
    error('iat_feat_based_init.m: Input "method" must be a string');
end

if (~(strcmp(transform,'homography')||strcmp(transform,'affine')||...
      strcmp(transform,'similarity')||strcmp(transform,'euclidean')||...
      strcmp(transform,'translation')))
  error('iat_feat_based_init.m: Input "transform" must be a valid string');
end

if (~(strcmp(method,'LS')||strcmp(method,'RANSAC')||...
      strcmp(method,'ext-LS')||strcmp(method,'ext-RANSAC')||...
      strcmp(method,'desc-LS')||strcmp(method,'desc-RANSAC')))
  
  error('iat_feat_based_init.m: Input "method" must be a valid string');
end

% Number of varargin inputs
vars=size(varargin);

%% Extract required inputs
if (strcmp(method,'LS')||strcmp(method,'RANSAC'))
    if vars(2)<2
        error('iat_feat_based_init.m: Image and/or template missing');
    end
    
    image=varargin{1};
    template=varargin{2};
    offset=3;
elseif (strcmp(method,'ext-RANSAC')||strcmp(method,'ext-LS'))
    if vars(2)<2
        error('iat_feat_based_init.m: Image and/or template points missing');
    end
    
    % Extract image and template points
    iPts=varargin{1};
    tPts=varargin{2};
    
    if ( (size(iPts,1)~=2) && (size(iPts,1)~=3) )
        error('iat_feat_based_init.m: image points structure must be 2xN or 3xN');
    end
    
    if ( (size(tPts,1)~=2) && (size(tPts,1)~=3) )
        error('iat_feat_based_init.m: template points structure must be 2xN or 3xN');
    end
    
    % Pad points with scales if necessary
    if (size(iPts,1)==2)
        iPts=iat_scale_pad(iPts);
    end
    if (size(tPts,1)==2)
        tPts=iat_scale_pad(tPts);
    end
    
    if ( size(iPts,2)~=size(tPts,2) )
        error('iat_feat_based_init.m: image and template points must have same cardinality');
    end
    
    offset=3;
elseif (strcmp(method,'desc-RANSAC')||strcmp(method,'desc-LS'))
    if vars(2)<2
        error('iat_feat_based_init.m: Image and/or template descriptors missing');
    end
    
    iDesc=varargin{1};
    tDesc=varargin{2};
    
    if vars(2)<4
        error('iat_feat_based_init.m: Image and/or template points missing');
    end
    
    iPts=varargin{3};
    tPts=varargin{4};
    
    if ( (size(iPts,1)~=2) && (size(iPts,1)~=3) )
        error('iat_feat_based_init.m: image points structure must be 2xN or 3xN');
    end
    
    if ( (size(tPts,1)~=2) && (size(tPts,1)~=3) )
        error('iat_feat_based_init.m: template points structure must be 2xN or 3xN');
    end
    
    if (size(iPts,1)==2) iPts=iat_scale_pad(iPts), end
    if (size(tPts,1)==2) tPts=iat_scale_pad(tPts), end
    
    offset=5;
end

%% Optional input parsing

p=inputParser;
if (strcmp(method,'RANSAC'))
    addOptional(p,'rMaxIter',100,@(x) (isscalar(x) && x>0));
    addOptional(p,'rMaxInvCount',10,@(x) (isscalar(x) && x>0));
    addOptional(p,'rOutFreeProb',.99,@(x) (isscalar(x) && x>0 && x<1));
    addOptional(p,'rErrorThr',1,@(x) (isscalar(x) && x>0));
    addOptional(p,'distRatio',.9,@(x) (isscalar(x) && x>0));
    addOptional(p,'descSize',64,@(x) (isscalar(x) && (x==32)||(x==64)||x==128));
    
    parse(p,varargin{3:end});
    
    ransacIter=p.Results.rMaxIter;
    ransacInvCount=p.Results.rMaxInvCount;
    ransacProb=p.Results.rOutFreeProb;
    ransacErrorThr=p.Results.rErrorThr;
    distRatio=p.Results.distRatio;
    sizeD=p.Results.descSize;
    
elseif (strcmp(method,'ext-RANSAC'))
    addOptional(p,'rMaxIter',100,@(x) (isscalar(x) && x>0));
    addOptional(p,'rMaxInvCount',10,@(x) (isscalar(x) && x>0));
    addOptional(p,'rOutFreeProb',.99,@(x) (isscalar(x) && x>0 && x<1));
    addOptional(p,'rErrorThr',1,@(x) (isscalar(x) && x>0));
    
    parse(p,varargin{3:end});
    
    ransacIter=p.Results.rMaxIter;
    ransacInvCount=p.Results.rMaxInvCount;
    ransacProb=p.Results.rOutFreeProb;
    ransacErrorThr=p.Results.rErrorThr;
    
elseif (strcmp(method,'desc-RANSAC'))
    addOptional(p,'rMaxIter',100,@(x) (isscalar(x) && x>0));
    addOptional(p,'rMaxInvCount',10,@(x) (isscalar(x) && x>0));
    addOptional(p,'rOutFreeProb',.99,@(x) (isscalar(x) && x>0 && x<1));
    addOptional(p,'rErrorThr',1,@(x) (isscalar(x) && x>0));
    addOptional(p,'distRatio',.9,@(x) (isscalar(x) && x>0));
    
    parse(p,varargin{5:end});
    
    ransacIter=p.Results.rMaxIter;
    ransacInvCount=p.Results.rMaxInvCount;
    ransacProb=p.Results.rOutFreeProb;
    ransacErrorThr=p.Results.rErrorThr;
    distRatio=p.Results.distRatio;
    
elseif (strcmp(method,'LS'))
    addOptional(p,'distRatio',.9,@(x) (isscalar(x) && x>0));
    addOptional(p,'descSize',64,@(x) (isscalar(x) && (x==32)||(x==64)||x==128));
    parse(p,varargin{3:end});
    distRatio=p.Results.distRatio;
    sizeD=p.Results.descSize;
    
elseif (strcmp(method,'desc-LS'))
    addOptional(p,'distRatio',.9,@(x) (isscalar(x) && x>0));
    parse(p,varargin{3:end});
    distRatio=p.Results.distRatio;
end


%% Apply the appropriate algorithm
if (strcmp(method,'RANSAC')||strcmp(method,'ext-RANSAC')||strcmp(method,'desc-RANSAC'))
    % RANSAC CASE
    
    if strcmp(method,'RANSAC')  % Extract SURF-descriptors

        [x1 x2]=iat_get_correspondences(image,template,sizeD,distRatio);
        x1=iat_scale_pad(x1);
        x2=iat_scale_pad(x2);

    elseif strcmp(method,'ext-RANSAC') % External correspondencies
    
        x1=iPts;
        x2=tPts;
    
    elseif strcmp(method,'desc-RANSAC') % External descriptors
        
        % Try to find an appropriate matching for every template descriptor
        % and find the valid ones
        mapping=iat_match_descriptors(tDesc',iDesc',distRatio);
        [v i]=find(mapping>0);
        
        % Points according to mapping (indexes)
        TempPts=i;
        ImPts=mapping(TempPts);
        
        x1=iPts(:,ImPts);
        x2=tPts(:,TempPts);
        
        if (size(x1,1)==2) x1=iat_scale_pad(x1); end
        if (size(x2,1)==2) x2=iat_scale_pad(x2); end
        
    end
    
    % Apply RANSAC algorithm
    [inliers H]=iat_ransac(x2,x1,ransacProb,transform,ransacIter,ransacInvCount,ransacErrorThr);

elseif (strcmp(method,'LS')||strcmp(method,'ext-LS')||strcmp(method,'desc-LS'))
    % LEAST SQUARES CASE
    
    if strcmp(method,'LS')
        % Size of SURF descriptors
        sizeD = 64;

        % Convert main image and template to a suitable form
        if max(image(:))<=1
            image = double(image)*255;
        end
        if max(template(:))<=1
            template = double(template)*255;
        end

        % Extract descriptors from main image and template
        [descI,locI]=iat_surf_desc(image,sizeD);
        [descT,locT]=iat_surf_desc(template,sizeD);

        % Initialize correspondence array
        correspondencies=[];

        % Descriptor matching
        mapping=iat_match_descriptors(descT',descI',distRatio);
        [v i]=find(mapping>0);
        TempPts=i;
        ImPts=mapping(TempPts);
        correspondencies=[locI(ImPts,1:2) locT(TempPts,1:2)];
        
        % Add scales to points
        x1=iat_scale_pad(correspondencies(:,1:2)');
        x2=iat_scale_pad(correspondencies(:,3:4)');
        
    elseif strcmp(method,'ext-LS')
        x1=iPts;
        x2=tPts;
    elseif strcmp(method,'desc-LS')
        % Try to find an appropriate matching for every template descriptor
        % and find the valid ones
        mapping=iat_match_descriptors(tDesc',iDesc',distRatio);
        [v i]=find(mapping>0);
        
        % Points according to mapping (indexes)
        TempPts=i;
        ImPts=mapping(TempPts);
        
        x1=iPts(:,ImPts);
        x2=tPts(:,TempPts);
    end
    
    % Find an estimation according to given transform
    switch transform
        case 'homography'
            H=iat_get_homography(x1,x2);
        case 'affine'
            H=iat_get_affine(x1,x2);
        case 'similarity'
            H=iat_get_similarity(x1,x2);
        case 'euclidean'
            H=iat_get_euclidean(x1,x2);
        case 'translation'
            H=iat_get_translation(x1,x2);
    end
end

% Output warp
warp=H;