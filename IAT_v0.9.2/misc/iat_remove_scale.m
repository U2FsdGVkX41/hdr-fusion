function [ npoints ] = iat_remove_scale( points )
% NPOINTS = IAT_REMOVE_SCALE( POINTS )
% IAT_REMOVE_SCALE removes the scale from homogeneous coordinates in
% array POINTS by normalizing with the third coordinate. 
% Points at infinity are left unchanged.
%
% -->Input
% POINTS:                    3xN array of N points 
%
% -->Output
% NPOINTS:                   3xN array with normalized homogeneous coordinates
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

% Check for valid data points array
if (size(points,1)~=3)
    error('iat_homogeneous_coords: Input array must be 3xN');
end

% Initialize output
samples=size(points,2);
npoints=points;

% Bring all data points to a common scale of 1
for i=1:samples
   if (points(3,i)==0)
       warning('iat_homogeneous_coords: points at infinity detected!');
   else
       npoints(:,i)=points(:,i)/points(3,i);
   end
end

end