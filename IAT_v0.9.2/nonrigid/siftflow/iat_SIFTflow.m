function [vx,vy,energylist]=iat_SIFTflow(im1,im2,par)
% [VX, VY,ENERGY]=iat_SIFTflow(SIM1,SIM2,PAR)
% iat_SIFTflow computes the flows from sift-image SIM1 to sift-image SIM2
% based on SIFTflow algorithm [1]. SIFTflow estimates dense correspondences
% from SIM1 to SIM2 in a coarse-to-fine manner(multi-resolution
% implementation) is used here. SIFT images are images whose pixels are
% SIFT descriptors [2]. The images must not necessarily be of the same size.
%
% Input arguments:
% SIM1:                       the source SIFT-image
% SIM2:                       the target SIFT-image
% PAR:                        a parameter struct with the following fields:
%
%      PAR.alpha:             the scale of the truncated L1-norm regularizer
%                             of flow discontinuity (default: 0.01)
%      PAR.d:                 the threshold of the truncated L1-norm
%                             (default: PAR.alpha*20)
%      PAR.gamma:             the scale of the regularization on flow
%                             magnitude (L1 norm) (default: 0.001)
%      PAR.nlevels:           the number of levels of the Gaussian pyramid
%                             (default:4)
%      PAR.topwsize:          the size of the matching window at the top
%                             (coarsest) level (default: 10)
%      PAR.nTopIterations:    the number of Belief-Propagation iterations at
%                             the top (coarsest) level(default: 100)
%      PAR.wsize:             the size of the matching window at lower
%                             levels (default: 3)
%      PAR.nIterations:       the number of BP iterations at lower levels
%                             (default: 40)
%
% Output arguments
% VX:                         the horizontal components of flows (displacements)
% VY:                         the vertical componenets of flows
% ENERGY                      a struct array with the minimum energy
%                             achieved per iteration per level.
%                             ENERGY(i).data is a vector with PAR.nIterations
%                             elements, when i<nlevels and with PAR.nTopIterations
%                             elements when i=nlevels.
%
%
% References:
% [1] C. Liu, J. Yuen, A. Torralba: SIFT Flow: Dense Correspondence across
% Scenes and its Applications, IEEE Trans. on PAMI, vol. 33, no. 5, 2011
% [2] D.G. Lowe, Distinctive image features from scale-invariant keypoints,
% Int. Journal on Computer Vision, vol. 60, no. 2, 2004
%
% -------------------
% Authors: Ce Liu, Georgios Evangelidis
% Copyright (C) 2013 Ce Liu, CSAIL MIT, celiu@mit.edu
% All rights reserved.
%
% This function is a minor modification of the original function written by
% Ce Liu. For any bugs, please contact <celiu@microsoft.com> or
% <georgios@iatool.net>
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default parameters
par0.alpha = 0.01;
par0.gamma = 0.001;
par0.nlevels = 4;
par0.wsize = 3;
par0.topwsize = 10;
par0.nIterations = 40;
par0.nTopIterations = 100;


if exist('par','var')
    if ~isstruct(par)
        error('iat_SIFTflow: parameter structure is not a matlab struct');
    end
    % merge default with given values
    params = iat_merge_param(par0, par);
else
    params = par0;
end

alpha = params.alpha;
gamma = params.gamma;
nlevels = params.nlevels;
wsize = params.wsize;
topwsize = params.topwsize;
nIterations = params.nIterations;
nTopIterations = params.nTopIterations;


if isfield(params,'d')
    d=params.d;
else
    d=alpha*20;
end

if ~isfloat(im1)
    im1=im2double(im1);
end
if ~isfloat(im2)
    im2=im2double(im2);
end

% build the Gaussian pyramid for the SIFT images
pyrd(1).im1=im1;
pyrd(1).im2=im2;

for i=2:nlevels
    pyrd(i).im1=imresize(imfilter(pyrd(i-1).im1,fspecial('gaussian',5,0.67),'same','replicate'),0.5,'bicubic');
    pyrd(i).im2=imresize(imfilter(pyrd(i-1).im2,fspecial('gaussian',5,0.67),'same','replicate'),0.5,'bicubic');
end

for i=1:nlevels
    [height,width,nchannels]=size(pyrd(i).im1);
    [height2,width2,nchannels]=size(pyrd(i).im2);
    [xx,yy]=meshgrid(1:width,1:height);
    pyrd(i).xx=round((xx-1)*(width2-1)/(width-1)+1-xx);
    pyrd(i).yy=round((yy-1)*(height2-1)/(height-1)+1-yy);
end

%nIterationArray=round(linspace(nIterations,nIterations*0.6,nlevels));
nIterationArray = nIterations*ones(1,nlevels);
if nlevels==1
    fprintf('SIFTflow is running in single-level mode\n');
else
    fprintf('SIFTflow is running in multi-level mode\n');
end
for i=nlevels:-1:1
    if nlevels>1
        fprintf('Level: %d....',i);
    end
    [height,width,nchannels]=size(pyrd(i).im1);
    [height2,width2,nchannels]=size(pyrd(i).im2);
    [xx,yy]=meshgrid(1:width,1:height);
    
    if i==nlevels % top level
        vx=pyrd(i).xx;
        vy=pyrd(i).yy;
        
        winSizeX=ones(height,width)*topwsize;
        winSizeY=ones(height,width)*topwsize;
    else % lower levels
        vx=round(pyrd(i).xx+imresize(vx-pyrd(i+1).xx,[height,width],'bicubic')*2);
        vy=round(pyrd(i).yy+imresize(vy-pyrd(i+1).yy,[height,width],'bicubic')*2);
        
        winSizeX=ones(height,width)*wsize;
        winSizeY=ones(height,width)*wsize;
    end
    
    %     if nchannels<=3
    %         Im1=im2feature(pyrd(i).im1);
    %         Im2=im2feature(pyrd(i).im2);
    %     else
    %         Im1=pyrd(i).im1;
    %         Im2=pyrd(i).im2;
    %     end
    
    Im1=pyrd(i).im1;
    Im2=pyrd(i).im2;
    
    
    % call SIFTflow MEX file
    if i==nlevels
        [flow,foo]=iat_SIFTflow_mex(Im1,Im2,[alpha,d,gamma*2^(i-1),nTopIterations,2,topwsize],vx,vy,winSizeX,winSizeY);
    else
        [flow,foo]=iat_SIFTflow_mex(Im1,Im2,[alpha,d,gamma*2^(i-1),nIterationArray(i),nlevels-i,wsize],vx,vy,winSizeX,winSizeY);
    end
    energylist(i).data=foo;
    vx=flow(:,:,1);
    vy=flow(:,:,2);
    if nlevels>1
        fprintf('done\n');
    end
end

