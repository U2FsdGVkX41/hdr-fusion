function new_warp = iat_warp_update(warp,delta_p,transform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW_WARP = IAT_WARP_UPDATE(WARP,DELTA_P,TRANSFORM)
% IAT_WARP_UPDATE updates the warp parameters by adding the correction values
% of DELTA_P to the current WARP.
%
% -->Input:
% WARP:                 the current warp transform,
% DELTA_P:              the current correction parameter vector,
% TRANSFORM:            the type of adopted transform, valid strings:
%                       {'translation','euclidean','affine','homography'}.
%
% -->Output:
% NEW_WARP:             the new (updated) warp transform
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

if strcmp(transform,'homography')
    delta_p=[delta_p; 0];
    new_warp=warp + reshape(delta_p, 3, 3);
    new_warp(3,3)=1;
elseif strcmp(transform,'affine')
    new_warp(1:2,:)=warp(1:2,:) + reshape(delta_p, 2, 3);
elseif strcmp(transform,'translation')
    new_warp =warp + delta_p;
elseif strcmp(transform, 'euclidean')

    theta = sign(warp(2,1))*acos(warp(1,1))+delta_p(1);

    tx = warp(1,3)+delta_p(2);
    ty = warp(2,3)+delta_p(3);
    new_warp = [cos(theta) -sin(theta) tx;...
                sin(theta) cos(theta) ty];
else
    disp('iat_warp_update: Unknown transform name');
    
end
