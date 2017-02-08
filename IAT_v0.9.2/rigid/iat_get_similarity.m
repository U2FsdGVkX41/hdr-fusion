function [ H ] = iat_get_similarity( pts1, pts2 )
% H = IAT_GET_SIMILARITY( POINTS1, POINTS2 )
%
% IAT_GET_SIMILARITY computes the similarity transform H of size 3x3 
% so that POINTS2 = H*POINTS1 in the least-square sense
%
% -->Input:
% POINTS1:              2xN or 3xN array of N initial points 
% POINTS2:              2xN or 3xN array of N transformed points 
%
% -->Output:
% H:                    The similarity transform obtained from 
%                       correspondences POINTS1<-->POINTS2
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

if (nargin~=2)
    error('iat_get_similarity.m: Wrong number of arguments');
end

samples=size(pts1,2);

% The structure of both point inputs must be of same size
if ~all(size(pts1)==size(pts2))
    error('iat_get_similarity: point structures must be of same size');
end

% homogeneous points
if size(pts1,1)==2
    pts1 = [pts1;ones(1,samples)];
    pts2 = [pts2;ones(1,samples)];
end

% Normalizing input data points
[xA, TA]=iat_normal_coords(pts1);
[xB, TB]=iat_normal_coords(pts2);

% Initializing necessary arrays to solve the transform in the least squares
% sense
A=zeros(2*samples,4);
xp=zeros(2*samples,1);

if (samples>1)
    % Forming matrix A (A*h=xp)
    % h = [c s tx ty]
    for k=1:samples
        A(2*k-1,:)=[xA(1,k) -xA(2,k) 1 0];
        A(2*k,:)  =[xA(2,k) xA(1,k) 0 1];
        
        xp(2*k-1)=xB(1,k);
        xp(2*k)  =xB(2,k);
    end
    
    % Solving the system via SVD and denormalize extracted transform
    [U, D, V]=svd(A);
    h=V*pinv(D)*U'*xp;    
    H=[h(1) -h(2) h(3);h(2) h(1) h(4);0 0 1];
    H=TB\H*TA;
else
    error('iat_get_similarity: two samples are at least required to obtain a similarity transform');
end

end