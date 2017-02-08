function [ H ] = iat_get_translation( ptsA, ptsB )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% H = IAT_GET_TRANSLATION( POINTS1, POINTS2 )
%
% IAT_GET_TRANSLATION computes the translation transform H so that 
% POINTS2 = H*POINTS1 in the least-L1 sense (median)
%
% -->Input:
% POINTS1:              2xN or 3xN array of N initial points 
% POINTS2:              2xN or 3xN array of N transformed points 
%
% -->Output:
% H:                    The translation transform obtained from 
%                       correspondences POINTS1<-->POINTS2. H is a 3x3
%                       array
%                              1    0    Tx
%                        H = [ 0    1    Ty ]
%                              0    0    1
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
    error('iat_get_translation.m: Wrong number of arguments');
end

% The structure of both point inputs must be of same size
if ~all(size(ptsA)==size(ptsB))
    error('iat_get_translation: point structures must be of same size');
end

dif=ptsB-ptsA;
mDifX=median(dif(1,:));
mDifY=median(dif(2,:));

H=eye(3);
H(1,3)=mDifX;
H(2,3)=mDifY;

end

