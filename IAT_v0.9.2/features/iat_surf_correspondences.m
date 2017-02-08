function [ Points1, Points2 ] = iat_surf_correspondences( image1, image2, angleRatio )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [POINTS1, POINTS2] = IAT_SURF_CORRESPONDENCES( IMAGE1, IMAGE2, ANGLE_RATIO)
% IAT_SURF_CORRESPONDENCES returns the coordinates of matched keypoints in 
% images IMAGE1 and IMAGE2. After extracting SURF descriptors from images,
% it obtains matches so that minimumArcCosine<ANGLE_RATIO*secondMinimumArcCosine
% (inverse cosine of normalized correlation between vectors is used)
%
% -->Input:
% IMAGE1:               First image
% IMAGE2:               Second image
% ANGLE_RATIO:          Angle ratio for matching descriptors (see the
%                       help file of the function IAT_MATCH_DESCRIPTORS)
%
% -->Output:
% POINTS1:              Points in IMAGE1
% POINTS2:              Correspondences of POINTS1 in IMAGE2
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Convert input images to a valid format
if max(image1(:))<=1
    image1 = double(image1)*255;
end
if max(image2(:))<=1
    image2 = double(image2)*255;
end

% Extract SURF descriptors
[desc1,loc1]=iat_surf(image1, 128);
[desc2,loc2]=iat_surf(image2, 128);

% Match descriptors
[~, ~, index1, index2] = iat_match_features(desc1,desc2,angleRatio);
% Use the mex implementation below for efficiency
%[~, ~, index1, index2] = iat_match_features_mex(desc1,desc2,angleRatio);


% Final points (x1<->x2)
Points1 = loc1(index1,1:2);
Points2 = loc2(index2,1:2);

end

