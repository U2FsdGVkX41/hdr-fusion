function [ H ] = iat_get_homography( pts1, pts2 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ H] = IAT_GET_HOMOGRAPHY( POINTS1, POINTS2 )
% IAT_GET_HOMOGRAPHY computes the homography H from points correspondences
% POINTS1 <--> POINTS2, so that POINTS2 = H*POINTS1 using the 
% Direct Linear Transform (DLT) algorithm
%
% -->Input:
% POINTS1:              2xN or 3xN array of N initial points 
% POINTS2:              2xN or 3xN array of N transformed points 
%
% -->Output:
% H:                    The extracted homography
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
    error('iat_get_homography: Wrong number of arguments');
end

% The structure of both point inputs must be of same size
if ~all(size(pts1)==size(pts2))
    error('iat_get_homography: point structures must be of same size');
end

% Number of input points
points=size(pts1,2);

% homogeneous points if needed
if size(pts1,1)==2
    pts1 = [pts1;ones(1,points)];
    pts2 = [pts2;ones(1,points)];
end


% Normalizing the input coordinates
[pts1, TA]=iat_normal_coords(pts1);
[pts2, TB]=iat_normal_coords(pts2);

% Assembling A matrix which consists of all the necessery equations to
% obtain the final homography
A=zeros(3*points,9);
for i=1:points
    % Declaring variables for convenience
    xVec=pts1(:,i)';
    x=pts2(1,i);
    y=pts2(2,i);
    z=pts2(3,i);
    
    % Calculating Ai matrix and putting it to the final A matrix
    A(3*i-2,:)=[zeros(1,3) -z*xVec y*xVec];
    A(3*i-1,:)=[z*xVec zeros(1,3) -x*xVec];
    A(3*i,:)  =[-y*xVec x*xVec zeros(1,3)];
end

% Calculating the SVD decomposition of the A matrix and extracting the
% homography as the final column of the V matrix
[~, ~, V]=svd(A,0);
H=reshape(V(:,9),3,3)';
H=TB\H*TA;

H=H/H(3,3);

end