function J = iat_warp_jacobian(nx, ny, warp, transform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JAC = IAT_WARP_JACOBIAN(NX, NY, WARP, TRANSFORM)
% IAT_WARP_JACOBIAN computes jacobian JAC of transform WARP with respect
% to parameters. In case of homography/euclidean transform, J depends on
% the parameter values, while in affine/translation case is totally invariant.
%
% -->Input:
% NX:                   the x-coordinate values of the horizontal side of
%                       ROI (i.e. [xmin:xmax]),
% NY:                   the y-coordinate values of vertical side of ROI
%                       (i.e. [ymin:ymax]),
% WARP:                 the warp transform (used only in homography and
%                       euclidean case),
% TRANSFORM:            the type of adopted transform. Valid string:
%                       {'affine''homography','translation','euclidean'}
%
% -->Output:
% JAC:                  The jacobian J. The size of J is 
%                       2*length(NY) x NOP*length(NX), where NOP the number
%                       of warp parameters
%                       homography: 8, 
%                       affine: 6, 
%                       euclidean: 3,
%                       translation: 2
%
% This function is used from ECC and LucasKanade algorithms. Use this 
% function carefully out of this framework
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

snx=length(nx);
sny=length(ny);

Jx=nx(ones(1,sny),:);
Jy=ny(ones(1,snx),:)';
J0=0*Jx;
J1=J0+1;


switch lower(transform)
    case 'homography'
        
        xy=[Jx(:)';Jy(:)';ones(1,snx*sny)];
        
        
        %3x3 matrix transformation
        A = warp;
        A = A./A(3,3);
        
        % new coordinates
        xy_prime = A * xy;
        
        % division due to homogeneous coordinates
        xy_prime(1,:) = xy_prime(1,:)./xy_prime(3,:);
        xy_prime(2,:) = xy_prime(2,:)./xy_prime(3,:);
        
        den = xy_prime(3,:)';
        
        Jx(:) = Jx(:) ./ den;
        Jy(:) = Jy(:) ./ den;
        J1(:) = J1(:) ./ den;
        
        Jxx_prime = Jx;
        Jxx_prime(:) = Jxx_prime(:) .* xy_prime(1,:)';
        Jyx_prime = Jy;
        Jyx_prime(:) = Jyx_prime(:) .* xy_prime(1,:)';
        
        Jxy_prime = Jx;
        Jxy_prime(:) = Jxy_prime(:) .* xy_prime(2,:)';
        Jyy_prime = Jy;
        Jyy_prime(:) = Jyy_prime(:) .* xy_prime(2,:)';
        
        
        J = [Jx, J0, -Jxx_prime, Jy, J0, - Jyx_prime, J1, J0;...
            J0, Jx, -Jxy_prime, J0, Jy, -Jyy_prime, J0, J1];
        
    case 'affine'
        
        % warp is not taken into account
        
        J = [Jx, J0, Jy, J0, J1, J0;...
            J0, Jx, J0, Jy, J0, J1];
        
    case 'translation'
        
        % warp is not taken into account
        J = [J1, J0;...
            J0, J1];
        
    case 'euclidean'
        
        mycos = warp(1,1);
        mysin = warp(2,1);
        
        Jx_prime = -mysin*Jx - mycos*Jy;
        Jy_prime =  mycos*Jx - mysin*Jy;
        
        J = [Jx_prime, J1, J0;...
            Jy_prime, J0, J1];
        
    otherwise
        error('iat_warp_jacobian: Unknown transform name!');
end




