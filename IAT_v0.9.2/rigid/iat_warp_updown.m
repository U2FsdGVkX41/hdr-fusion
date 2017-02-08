function warp = iat_warp_updown(warp_in, transform, flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WARPNEW = IAT_WARP_UPDOWN(WARPOLD, TRANSFORM, FLAG)
% iat_warp_updown appropriately modifies the WARPOLD and returns WARPNEW
% in order to apply the latter in the next level of pyramid. If FLAG is
% equal to 1, the function makes WARPNEW appropriate for the down level of
% higher resolution. If FLAG is equal to 0, the function makes WARPNEW
% appropriate for the up level of lower resolution.
%
% -->Input:
% WARPOLD:              the current warp transform
% TRANSFORM:            the type of adopted transform. Valid strings:
%                       {'translation','euclidean','affine','homography'}.
% FLAG:                 The flag which defines the 'next' level. 1 means
%                       that the the next level is a higher resolution
%                       level, while 0 means that it is a lower resolution
%                       level.
%
% -->Output:
% WARPNEW:              the next-level warp transform
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

if ~iat_is_transform(transform)
    error('iat_warp_updown: non-valid transform name');
end

if strcmpi(transform,'similarity');
    error(['iat_warp_updown: current version of function does not support' ...
        ' similarity. Choose affine or euclidean instead']);
end

warp=warp_in;
if flag==1
    if strcmp(transform,'homography')
        warp(7:8)=warp(7:8)*2;
        warp(3)=warp(3)/2;
        warp(6)=warp(6)/2;
    end
    
    if strcmp(transform,'affine')
        warp(1:2,3) = warp(1:2,3)*2;
        
    end
    
    if strcmp(transform,'translation')
        warp = warp*2;
    end
    
    if strcmp(transform,'euclidean')
        warp(1:2,3) = warp(1:2,3)*2;
    end
    
end

if flag==0
    if strcmp(transform,'homography')
        warp(7:8)=warp(7:8)/2;
        warp(3)=warp(3)*2;
        warp(6)=warp(6)*2;
    end
    
    if strcmp(transform,'affine')
        warp(1:2,3) = warp(1:2,3)/2;
    end
    
    if strcmp(transform,'euclidean')
        warp(1:2,3) = warp(1:2,3)/2;
    end
    
    if strcmp(transform,'translation')
        warp = warp/2;
    end
    
end