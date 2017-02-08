clc;
close all;
clear all;

%% Load and show images
% Read initial image <-> template pair
img=imread('G:\HDRImageSynthesis\PatchBasedHDR_MATLAB_v1.1\Scenes\lady\Img1.tif');
tmp=imread('G:\HDRImageSynthesis\PatchBasedHDR_MATLAB_v1.1\Scenes\lady\Img2.tif');


transform = 'homography'; % the ground transformation is LensUndistortion+Homography (images are real)

% Plot images
figure;imshow(img);title('Image', 'Fontsize', 14);
figure;imshow(tmp);title('Template', 'Fontsize', 14);

%% Compute SURF descriptors and show correspondences
% Extract SURF descriptors
[d1, l1]=iat_surf(img,128);
[d2, l2]=iat_surf(tmp,128);

% An initial match of the descriptors
[map, match, imgInd, tmpInd]=iat_match_features(d1,d2,.8);
iat_plot_correspondences(img,tmp,l1(imgInd,1:2)',l2(tmpInd,1:2)');
title('RANSAC-filtered correspondences', 'Fontsize', 14);


% because of the chess content, the correspondences are geometrically
% consistent and feature-based alignment is meaningless


%% Area-based alignment with no initialization
par = [];
par.transform = transform;
par.levels = 1;
par.iterations = 30;

% run ECC (you can try here Lucas-Kanade method also)
[ECCWarp]=iat_ecc(img,tmp,par); 

% Compute the warped image and visualize the error
[wimageECC, supportECC] = iat_inverse_warping(img, ECCWarp, par.transform, 1:size(tmp,2),1:size(tmp,1));

%plot the warped image
figure;imshow(uint8(wimageECC)); title('Warped image by ECC', 'Fontsize', 14);

% visualize the error
[imdiffECC, grayerrorECC] = iat_error2gray(tmp,wimageECC,supportECC);
figure;imshow(grayerrorECC); title('Error of ECC alignment', 'Fontsize', 14);

% Mosaic Image
C=iat_mosaic(tmp,img,ECCWarp);
figure;imshow(uint8(C));title('Mosaic Image using ECC', 'Fontsize', 14);