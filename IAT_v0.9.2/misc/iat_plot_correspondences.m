function iat_plot_correspondences( img1, img2, pts1, pts2 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IAT_PLOT_CORRESPONDENCES(IMAGE1, IMAGE2, POINTS1, POINTS2 )
% IAT_PLOT_CORRESPONDENCES draws lines between correspoding points of
% IMAGE1 and IMAGE2, i.e. POINTS1 and POINTS2. Drawing takes place in an
% image defined by the horizontal concatenation of IMAGE1 and
% IMAGE2.
%
% -->Input:
% IMAGE1:                 The first image
% IMAGE2:                 The second image
% POINTS1:                The points of IMAGE1 (2xN array)
% POINTS2:                The points of IMAGE2 that corresponds to POINTS1
%                         (2xN array)
%
% -->Output:
% NONE. The function just exports the appropriate plot.
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

% Extract image sizes
chA=size(img1);
chB=size(img2);

if numel(chA)==2
    chA(3) = size(img1,3);
end

if numel(chB)==2
    chB(3) = size(img2,3);
end

if chA(3)==1 && chB(3)==3
    img1 = img1(:,:,[1 1 1]);
    chA(3)=3;
end

if chA(3)==3 && chB(3)==1
    img2 = img2(:,:,[1 1 1]);
    chB(3)=3;
end


%Color format check
if (chA(3)~=chB(3))
    error('iat_plot_correspondences: Images must be in common color format');
end

% Point structures check
if (size(pts1,1)~=2) || (size(pts2,1)~=2)
    error('iat_plot_correspondences:point structures must be of size 2xN');
elseif (size(pts1,2)~=size(pts2,2))
    error('iat_plot_correspondences: point structures must have the same size');
end

% Number of points
points=size(pts1,2);

% Number of padding rows
padding=abs(chA(1)-chB(1));

% Pad image with smallest height with zeros
if (chA(1)>chB(1))
    if (chA(3)==1)
        padArray=zeros(padding,chB(2));
    elseif (chA(3)==3)
        padArray=zeros(padding,chB(2),3);
    else
        error('iat_plot_correspondences: Unknown color format in first image');
    end
    img2=[img2;padArray];
elseif (chB(1)>chA(1))
    if (chB(3)==1)
        padArray=zeros(padding,chA(2));
    elseif (chB(3)==3)
        padArray=zeros(padding,chA(2),3);
    else
        error('iat_plot_correspondences: Unknown color format in second image');
    end
    img1=[img1;padArray];
end

% Concatenate images horizontally
imgC=[img1 img2];

% Plot concatenated images
figure;imshow(imgC);

% Translate points of second image in order to correspond in the
% concatenated image
pts2(1,:)=pts2(1,:)+size(img1,2);

hold on;

% Plot keypoints
plot(pts1(1,:),pts1(2,:),'rx','MarkerSize',12);
plot(pts2(1,:),pts2(2,:),'rx','MarkerSize',12);

% Plot correspondence lines
for i=1:points
    line([pts1(1,i) pts2(1,i)],[pts1(2,i) pts2(2,i)],'LineWidth',0.3,'Color',[.1 1 .1]);
end

hold off;
%title('Point Correspondences', 'Fontsize', 14);
end

