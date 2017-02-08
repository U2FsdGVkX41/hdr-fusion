function [ Inliers, Model ] = iat_ransac2( pts1, pts2, T, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ INLIERS, MODEL ] = IAT_RANSAC( PTS1, PTS2, TRANSFORM)
% IAT_RANSAC implements the ransac algorithm adapted to estimate the
% optimum TRANSFORM from point correspondences PTS1<->PTS2 and 
% returns the largest inlier set of correspondences and the respective
% transformation. Note that it uses a symmetric error function to choose 
% the inliers.
%
% -->Input:
% PTS1:                 A 3xN array that holds the initial set of
%                       homogeneous coordinates
% PTS2:                 A 3xN array that holds the final set of homogeneous
%                       coordinates (in correspondence with PTS1)
% TRANSFORM:            The type of the transform. It must be one
%                       of the following strings: {'homography','affine',
%                       'similarity','euclidean','translation'}
%
% -->Optional Paramaters:
% [ INLIERS, MODEL ] = IAT_RANSAC( PTS1, PTS2, TRANSFORM, 'PARAM1', PARAM1VALUE,...)
% The user can define his own parameter values insted of default values.
% These parameters are:
%
% 'maxIter':            The maximum number of iterations allowed. You may
%                       need more iterations if the result is not
%                       sufficient (default: 100).
%
% 'prob':               The desired probability to pick a MSS(Minimum
%                       Sample Set) with no outliers (default: 0.99).
%
% 'tol':                The tolerance for the symmetric error that 
%                       determines whether a data point is an inlier or 
%                       an outlier. Since the coordinates of points are
%                       normalized before the RANSAC core steps, the error
%                       tolerance can be lower than one (default: 1).
%
% 'maxInvalidCount':    The maximum number of degenerated set picks
%                       (default: 100)
%
%
% -->Output:
% INLIERS:              The indices of points in PTS1 (or PTS2) that belong
%                       to the largest inlier set
% MODEL:                The RANSAC-optimum estimated transform. It 
%                       corresponds to the output of the function 
%                       iat_get_<TRANSFORM> when using the inlier set only.
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


% input parser
par = inputParser;
addOptional(par,'maxIter',100);
addOptional(par,'prob',0.99);
addOptional(par,'tol',0.5);
addOptional(par,'maxInvalidCount',100);
parse(par,varargin{:});

%rename struct values   
maxIterations = par.Results.maxIter;
P = par.Results.prob;
errorThr = par.Results.tol;
maxInvalidCount = par.Results.maxInvalidCount;


% structure of arrays pts1 and pts2
n1=size(pts1);
n2=size(pts2);


% Parameters validity checking
if (~isscalar(maxIterations))
    error('iat_ransac.m: Parameter "maxIter" must be scalar');

elseif (maxIterations<1)
    error('iat_ransac.m: Parameter "maxIter" must be greater than one');
    
elseif (~isscalar(P))
    error('iat_ransac.m: Parameter "prob" must be scalar');
    
elseif (P<=0 || P>=1)
    error('iat_ransac.m: Parameter "prob" must be in the interval (0,1)');
    
elseif (~isscalar(errorThr))
    error('iat_ransac.m: Parameter "tol" must be scalar');
    
elseif (errorThr<0)
    error('iat_ransac.m: Parameter "maxIter" must be greater than zero');
    
elseif (~isscalar(maxInvalidCount))
    error('iat_ransac.m: Parameter "maxInvalidCount" must be scalar');
    
elseif (maxInvalidCount<1);
    error('iat_ransac.m: Parameter "maxInvalidCount" must be greater than one');
    
end

% Check transform validity
if ~(iat_is_transform(T))
    error('iat_ransac: Input parameter "Transform" does not represent a valid transform');
end

% Check for structure validity of arrays pts1 and pts2
if ~all(n1==n2)
    error('iat_ransac: Structures pts1 and pts2 must be of same size');
end

% Check for scale component in correspondences
if (n1~=3)
    error('iat_ransac: Structures pts1 and pts2 must be 3xN');
end

% Normalise the input points such that their centroid is the origin (0,0)
% and their mean distance from the centroid is sqrt(2). If you avoid this
% step, you should increase the order of error tolerance 
[pts1, T1]=iat_normal_coords(pts1);
[pts2, T2]=iat_normal_coords(pts2);

% Number of given samples as input in the algorithm
samples=n1(2);

% Minimum number of correspondences needed to fully obtain a model. Depends
% on the input transform
if (strcmp(T,'homography'))
    minS=4;
elseif (strcmp(T,'affine'))
    minS=3;
elseif (strcmp(T,'euclidean')||strcmp(T,'similarity'))
    minS=2;
elseif (strcmp(T,'translation'))
    minS=1;
end

% Required samples check
if (samples<minS)
    error('iat_ransac: %d samples are at least required',minS);
end

% Output arguments:
% 1) The estimated model
% 2) The inliers to this model
Model=[];
Inliers=[];


% Auxiliary variables
bestConsensus=0;
iterations=Inf;

% Iteration counter
i=0;

% A flag that shows that the algorithm exceeded the maximum number of
% iterations
failFlag=0;

while (i<iterations)
    % Check for maximum iterations violation
    if (i>maxIterations)
        warning('iat_ransac: Maximum number of iterations reached. Exiting algorithm...');
        failFlag=1;
        break
    end
    
    % Degenerated configurations in a row count initialization
    iCount=0;
    breakFlag=1;
    
    while (iCount<maxInvalidCount && breakFlag)
        % Generate an MSS (Minimum Sample Set) - Indexes
        perm=randperm(samples);
        initSet=perm(1:minS);
        
        % MSS correspondences
        xA=pts1(:,initSet);
        xB=pts2(:,initSet);
        
        % Check for degeneracy
        if (minS>2)
            breakFlag = iat_collinearity_check(xA) || iat_collinearity_check(xB);
        end
        iCount=iCount+1;
    end
    
    % Compute model from MSS correspondences
    
    if (strcmp(T,'homography'))
        tmpH=iat_get_homography(xA, xB);
    elseif (strcmp(T,'affine'))
        tmpH=iat_get_affine(xA, xB);
    elseif (strcmp(T,'euclidean'))
        tmpH=iat_get_euclidean(xA, xB);
    elseif (strcmp(T,'similarity'))
        tmpH=iat_get_similarity(xA,xB);
    elseif (strcmp(T,'translation'))
        tmpH=iat_get_translation(xA,xB);
    end
    
    % Temporary consensus set. The inliers are marked with 1's
    tmpConSet=zeros(1,samples);
    
    % Measuring distances of datapoints from estimated model
    for j=minS+1:samples
        
        % j-th pair of correspondences
        initX  = iat_remove_scale(pts1(:,perm(j)));
        finalX = iat_remove_scale(pts2(:,perm(j)));
        
        % Scale normalization before error measurement
        Hx   = iat_remove_scale(tmpH*initX);

        iHxp = iat_remove_scale(tmpH\finalX);
        
        % Calculate error for current point (Symmetric Transfer Error)
        dF = sum((finalX-Hx).^2);
        dB = sum((initX-iHxp).^2);
        d  = dF+dB;
        
        % Determine whether data point is an inlier or an outlier
        if (d<errorThr)
            tmpConSet(perm(j))=1;
        end
        
    end
    
    % Determine whether the current model is better than the previous one.
    % This is done by checking the cardinality of the new consensus set
    tmpConsensus=nnz(tmpConSet);
    if (tmpConsensus>bestConsensus)
        
        % In case of a better estimation, update the consensus of the best
        % set, keep track of the new inliers and the new candidate model
        bestConsensus=tmpConsensus;
        Inliers=find(tmpConSet==1);
        Model=tmpH;

        % check division with -inf and 0??? 
        % Estimate the new number of necessary iterations
        denom=log(1-(tmpConsensus/samples)^minS);
        if (denom==0)
            iterations=maxIterations;
        elseif (denom==-Inf)
            iterations=0;
        else
            iterations=ceil(log(1-P)/denom);
        end
        %         iterations=ceil(log(1-P)/log(1-(tmpConsensus/samples)^minS));
        
    end
    
    % Iteration counter update
    i=i+1;
    
end

if (length(Inliers)>minS-1)
    % Estimate final model using only the found inliers and denormalize model
    if (strcmp(T,'homography'))
        % Homography Case
        Model=T2\(iat_get_homography(pts1(:,Inliers),pts2(:,Inliers)))*T1;
        
    elseif(strcmp(T,'affine'))
        % Affine Case
        Model=T2\(iat_get_affine(pts1(:,Inliers),pts2(:,Inliers)))*T1;
        
    elseif(strcmp(T,'similarity'))
        % Similarity Case
        Model=T2\(iat_get_similarity(pts1(:,Inliers),pts2(:,Inliers)))*T1;
        
    elseif(strcmp(T,'euclidean'))
        % Euclidean Case
        Model=T2\iat_get_euclidean(pts1(:,Inliers),pts2(:,Inliers))*T1;
        
    elseif(strcmp(T,'translation'))
        % Translation Case
        Model=T2\(iat_get_translation(pts1(:,Inliers),pts2(:,Inliers)))*T1;
    else
        if ~failFlag
            fprintf(1,'iat_ransac: Not enough inliers were found. Try another configuration\n');
            Model=eye(3);
        end
    end
    
    Model=Model/Model(3,3);
    Model = double(Model);
    
    % Executed iterations
    fprintf(1,'RANSAC: Done in %d Iterations\n',i);
    
end