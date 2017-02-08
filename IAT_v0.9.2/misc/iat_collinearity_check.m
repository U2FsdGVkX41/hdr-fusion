function [ check ] = iat_collinearity_check( pts )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ CHECK ] = IAT_COLLINEARITY_CHECK( POINTS )
% IAT_COLLINEARITY_CHECK checks whether a set of up to four two-dimensional 
% points contains three collinear points; homogeneous coordinates are valid
%
% -->Input
% POINTS:               2x3 or 2x4 array of points. When homogeneous
%                       coordinates are used, the arrays have size 3x3 or 
%                       3x4 respectively
%
% -->Output
% CHECK:               A logical variable showing if collinearity is detected
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

% Initially, the collinearity check is false
check=0;

if size(pts,1)==2
    pts=iat_homogeneous_coords(pts);
end

if size(pts,1)~=3
    error('iat_collinearity_check: an array of 2D points should be given');
end

% Number of samples
samples=size(pts,2);

if (samples==4)
    K=4;
    
    A=cell(4,1);
    
    A{1}=[pts(:,1)';pts(:,2)';pts(:,3)'];
    A{2}=[pts(:,1)';pts(:,2)';pts(:,4)'];
    A{3}=[pts(:,1)';pts(:,3)';pts(:,4)'];
    A{4}=[pts(:,2)';pts(:,3)';pts(:,4)'];
    
elseif(samples==3)
    K=1;
    
    A=cell(1);
    
    A{1}=[pts(:,1)';pts(:,2)';pts(:,3)'];
elseif (samples==2);
    K=0;
else
    error('iat_collinearity_check: Number of points must be 2, 3 or 4');
end

% Compute the determinants of matrices A and warn if a collinearity
% exists
for i=1:K
    if abs(det(A{i}))<1e-16
        check=1;
        warning('Collinearity detected');
        break;
    end
end

check = logical(check);

