function [ hpoints ] = iat_homogeneous_coords(points, dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HPOINTS = IAT_HOMOGENEOUS_COORDS(POINTS, DIM)
% IAT_HOMOGENEOUS_COORDS converts DIM-dimensional into homogeneous 
% coordinates by adding an extra coordinate at each point (column) of POINTS.
% To avoid rescaling all coordinates, the extra coordinate is equal to one. 
% DIM should be given when it comes to more than two dimensions.
%
% -->Input:
% POINTS:               DIMxN array with N DIM-dimensional points
% DIM:                  The optional dimensionality of points (default: 2)
%
% -->Output:
% HPOINTS:              (DIM+1)xN array that contains the points of POINTS with
%                       homogeneous coordinates.
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

if nargin==1
    dim = 2;
end

% Check data points array for valid dimensionality
if size(points,1)~=dim
    error(['iat_homogeneous_coords: the points should be ' num2str(dim) '-dimensional']);
end

% Number of points;
N=size(points,2);
hpoints=[points;ones(1,N)]; % homogeneous coordinates


end