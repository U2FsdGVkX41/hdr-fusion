function G = iat_image_jacobian(gx, gy, jac, nop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G = IAT_IMAGE_JACOBIAN(IX, IY, JAC, NOP)
% IAT_IMAGE_JACOBIAN computes the jacobian G of the warped image with 
% respect to transformation parameters. 
% This matrix depends on the gradient of the warped image, as 
% well as on the jacobian JAC of the warp transform wrt parameters. 
%
% For a detailed definition of matrix G, see [1]
%
% -->Input:
% IX:                   The warped image gradient in x (horizontal) direction
% IY:                   The warped image gradient in y (vertical) direction
% JAC:                  The jacobian matrix J of the warp transform with 
%                       respect to parameters. The size of JAC is MxN, 
%                       where N = NOP*size(IX,2) and M = size(IX,1)
% NOP:                  The number of warp parameters:
%                       homography: 8, 
%                       affine: 6, 
%                       euclidean: 3,
%                       translation: 2
%
% -->Output:
% G:                    The jacobian matrix G.
% 
% This function is used by ECC and Lucas-Kanade algorithms. Use this 
% function carefully out of this framework
%
% References
% [1] "G.D.Evangelidis, E.Z.Psarakis, Parametric Image Alignment using 
% Enhanced Correlation Coefficient.IEEE Trans. on PAMI, vol.30, no.10, 2008"
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

[h,w]=size(gx);

if nargin<4
    error('iat_image_jacobian: not enough input arguments');
end


if (size(gx,2)*nop~=size(jac,2)) || (size(gx,2)*nop~=size(jac,2))
        error('iat_image_jacobian: incosistent size of wapred gradients and warp-jacobian');
end


gx=repmat(gx,1,nop);
gy=repmat(gy,1,nop);

G=gx.*jac(1:h,:)+gy.*jac(h+1:end,:);
G=reshape(G,h*w,nop);