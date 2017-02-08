function iat_warping_area(image, warp, nx, ny)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IAT_WARPING_AREA(IMAGE, WARP, NX, NY)
% IAT_WARPING_AREA marks the image area that should be backward warped.
% This area is the map of [NY,NX] grid through matrix WARP.
%
% -->Input:
% IMAGE:                the image to be warped,
% WARP:                 the warp transform,
% NX:                   the x-coordinate values of horizontal side of ROI
%                       (i.e. [xmin:xmax]),
% NY:                   the y-coordinate values of vertical side of ROI
%                       (i.e. [ymin:ymax]),
%
% -->Output:
% NONE
%
% -------------------
% Authors: Georgios Evangelidis
% Copyright (C) 2013 Georgios Evangelidis
% All rights reserved.
%
% For any bugs, please contact <georgios.evangelidis@inria.fr> or
% <anatolitis@ceid.upatras.gr>
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if size(warp,2)==1 && size(warp,1)==2
    warp = [eye(2) warp];
    warp = [warp; 0 0 1];
end

if size(warp,2)==3 && size(warp,1)==2
    warp = [warp; 0 0 1];
end



x1 = [nx(1) ny(1) 1]';
x2 = [nx(1) ny(end) 1]';
x3 = [nx(end) ny(end) 1]';
x4 = [nx(end) ny(1) 1]';
xy = [x1 x2 x3 x4];



%3x3 matrix transformation
A = warp;

A = A./A(3,3);


% new coordinates
xy_prime = A * xy;


xy_prime(1,:) = xy_prime(1,:)./xy_prime(3,:);
xy_prime(2,:) = xy_prime(2,:)./xy_prime(3,:);


% Ignore third row
X = xy_prime(1:2,:);

lw = 2;
figure;imshow(uint8(image));hold on;
line([X(1,1) X(1,2)], [X(2,1) X(2,2)], 'Linewidth', lw, 'color', [1 0 1]);
line([X(1,2) X(1,3)], [X(2,2) X(2,3)], 'Linewidth', lw, 'color', [1 0 1]);
line([X(1,3) X(1,4)], [X(2,3) X(2,4)], 'Linewidth', lw, 'color', [1 0 1]);
line([X(1,4) X(1,1)], [X(2,4) X(2,1)], 'Linewidth', lw, 'color', [1 0 1]);
hold off
