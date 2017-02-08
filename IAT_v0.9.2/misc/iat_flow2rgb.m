function rgb = iat_flow2rgb( fx, fy, dispRange )
% RGBIM = IAT_FLOW2RGB( FX, FY, DISPRANGE )
%
% IAT_FLOW2RGB creates a color image RGBIM from the dense displacements (flows) 
% FX and FY. The angle of each displacement vector is translated to a color
% through a specific colomap. DISPRANGE is an optional parameter that defines 
% the allowable range of DISPLACEMENT (e.g. [min_flow, max_flow]). 
% 
% -->Input
% FX:           Array with horizontal displacements
% FY:           Array with vertical dispalcements
% DISPRANGE:    The allowable range of displacements (default: [-1e9, 1e9])
% 
% -->Output
% 
% RGBIM:        The RGB image that visualizes the displacements
% 
% NOTE:   This function was built based on Deqing Sun's flow visualization 
% function, while it uses exactly the same colormap with Deqing's function
% (see below the function makeColorwheel())
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



if ~all(size(fx)==size(fy))
    error('iat_flow2rgb:flow arrays must have the same size');
end
if (length(size(fx))~=2)
    error('iat_flow2rgb:flows must be two dimensional arrays');
end
if nargin==3
    if dispRange(2)<=dispRange(1)
        error('iat_flow2rgb: non valid flow range');
    end
end


if nargin<3
    fr = [-1e9 1e9];
else
    fr = dispRange;
end


outOfRange = fx>fr(2) | fy>fr(2) | fx<fr(1) | fy <fr(1) | isnan(fx) | isnan(fy);
fx(outOfRange) = 0;
fy(outOfRange) = 0;


% flow normalization
mag = sqrt(fx.^2+fy.^2);
maxmag = max(mag(:));

fx = fx/maxmag;
fy = fy/maxmag;

% Create palette of colors
%palette = iat_flow_palette();
palette = makeColorwheel();

colors  = size(palette,1);

mag = sqrt(fx.^2+fy.^2);
angle  = atan2(-fy,-fx)/pi;

% makes angles a color index map
angle = (angle+1)/2; % angles in [0 1]
angle = angle*(colors-1)+1;

k0=floor(angle);
k1=k0+1;
k1(k1==colors+1)=1;

f=angle-k0;

[row,col] = size(fx);

red = zeros(row,col);
gre = red;
blu = red;
idx = mag<=1;


red(:) = ((1-f(:)).*palette(k0,1)+f(:).*palette(k1,1))/255;
red(idx) = 1-mag(idx).*(1-red(idx));
red(~idx) = red(~idx).*0.75;

gre(:) = ((1-f(:)).*palette(k0,2)+f(:).*palette(k1,2))/255;
gre(idx) = 1-mag(idx).*(1-gre(idx));
gre(~idx) = gre(~idx).*0.75;

blu(:) = ((1-f(:)).*palette(k0,3)+f(:).*palette(k1,3))/255;
blu(idx) = 1-mag(idx).*(1-blu(idx));
blu(~idx) = blu(~idx).*0.75;


rgb = cat(3,red,gre,blu);
rgb(repmat(outOfRange,[1,1,3]))=0;
rgb = uint8(rgb*255);
end


function palette = iat_flow_palette()
% IAT color palette 

% Create a palette of fully saturated colors
palette(:,1)=0:1/54:1;
palette(:,2)=1;
palette(:,3)=1;

% Convert HSV representation to RGB
palette=floor(hsv2rgb(palette).*255);

end

function colorwheel = makeColorwheel()
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-10-31 21:20:30 (Wed, 31 Oct 2006) $

% Copyright 2007, Deqing Sun.

%   color encoding scheme

%   adapted from the color circle idea described at
%   http://members.shaw.ca/quadibloc/other/colint.htm


RY = 15;
YG = 6;
GC = 4;
CB = 11;
BM = 13;
MR = 6;

ncols = RY + YG + GC + CB + BM + MR;

colorwheel = zeros(ncols, 3); % r g b

col = 0;
%RY
colorwheel(1:RY, 1) = 255;
colorwheel(1:RY, 2) = floor(255*(0:RY-1)/RY)';
col = col+RY;

%YG
colorwheel(col+(1:YG), 1) = 255 - floor(255*(0:YG-1)/YG)';
colorwheel(col+(1:YG), 2) = 255;
col = col+YG;

%GC
colorwheel(col+(1:GC), 2) = 255;
colorwheel(col+(1:GC), 3) = floor(255*(0:GC-1)/GC)';
col = col+GC;

%CB
colorwheel(col+(1:CB), 2) = 255 - floor(255*(0:CB-1)/CB)';
colorwheel(col+(1:CB), 3) = 255;
col = col+CB;

%BM
colorwheel(col+(1:BM), 3) = 255;
colorwheel(col+(1:BM), 1) = floor(255*(0:BM-1)/BM)';
col = col+BM;

%MR
colorwheel(col+(1:MR), 3) = 255 - floor(255*(0:MR-1)/MR)';
colorwheel(col+(1:MR), 1) = 255;
end
