function [ nPts, T ] = iat_normal_coords( pts )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ NPOINTS, T ] = IAT_NORMAL_COORDS( POINTS )
% IAT_NORMAL_COORDS implements the normalization of the coordinates of point
% set POINTS, such that their centroid is at the origin (0,0) and 
% their mean distance from origin is sqrt(2). It is an essential 
% preprocessing step before any transform estimation for better 
% stability/accuracy
%
% -->Input:
% POINTS:                The input set of points in homogeneous coordinates 
%
% -->Output:
% NPOINTS:               The normalized set of points
% T:                     The normalization transform, so that 
%                        NPOINTS = T*POINTS
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

% Scale-normalised coordinates
pts(1,:)=pts(1,:)./pts(3,:);
pts(2,:)=pts(2,:)./pts(3,:);
pts(3,:)=1;

% Calculate the centroid of the coordinates and subtract it to shift them
% to the origin (0,0)
center=mean(pts(1:2,:),2);
tmpPts=pts;
tmpPts(1,:)=tmpPts(1,:)-center(1);
tmpPts(2,:)=tmpPts(2,:)-center(2);

% Calculate the mean distance from their new centroid
distances=sqrt(tmpPts(1,:).^2+tmpPts(2,:).^2);
mDist=mean(distances(:));

% The new scaling factor
s=sqrt(2)/mDist;

% The normalization process as a transformation matrix
T=[s 0 -s*center(1);0 s -s*center(2);0 0 1];
nPts=T*pts;

end