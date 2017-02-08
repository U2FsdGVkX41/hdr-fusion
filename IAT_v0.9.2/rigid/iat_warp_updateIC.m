function new_warp = iat_warp_updateIC(warp,delta_p,transform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW_WARP = iat_warp_updateIC(WARP,DELTA_P,TRANSFORM)
% iat_warp_updateIC composes WARP with the inverse of DELTA_P to produce the
% new warp NEW_WARP (used by inverese-compositional alignment schemes)
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
% Authors: Georgios Evangelidis,
% Copyright (C) 2013 Georgios Evangelidis
% All rights reserved.
% 
% For any bugs, please contact <georgios.evangelidis@inria.fr> or 
% 
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% delta_p is always the correction of the indentity warp applied to 
% the template. The template-warp is inverted and composed with 
% the image-warp.

if strcmp(transform,'homography')
    delta_p=[delta_p; 0];
    
    idelta = inv(eye(3)+reshape(delta_p,3,3));
    idelta = idelta./idelta(3,3);
    
    new_warp=idelta*warp;
    new_warp = new_warp/new_warp(3,3);

elseif strcmp(transform,'affine')
 
    
    delta_p = reshape(delta_p, 2, 3);
    delta_p(1,1) = delta_p(1,1)+1;
    delta_p(2,2) = delta_p(2,2)+1;
    
    idelta = [delta_p; 0 0 1];	
    
    new_warp = idelta\[warp(1:2,:); 0 0 1];
    new_warp = new_warp(1:2,:);
    
elseif strcmp(transform,'translation')
    new_warp = warp - delta_p;
    
elseif strcmp(transform, 'euclidean')

    dtheta = delta_p(1);
    
    idelta = [cos(dtheta) -sin(dtheta) delta_p(2);...
        sin(dtheta) cos(dtheta) delta_p(3);...
        0 0 1];
    
    new_warp  = idelta\[warp; 0 0 1];
    new_warp = new_warp(1:2,:);
    
else
    disp('iat_warp_updateIC: Unknown transform name');
    
end

