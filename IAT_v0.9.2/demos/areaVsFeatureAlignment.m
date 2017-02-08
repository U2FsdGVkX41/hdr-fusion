clc;
clear all;
close all;

% comment out the figures you do not want to show

%% Read initial image <-> template pair
img=imread('G:\HDRImageSynthesis\PatchBasedHDR_MATLAB_v1.1\Scenes\lady\Img3.tif');
tmp=imread('G:\HDRImageSynthesis\PatchBasedHDR_MATLAB_v1.1\Scenes\lady\Img2.tif');

transform = 'homography';

% Plot both of these images
figure;imshow(img);title('Image','Fontsize',14);
figure;imshow(tmp);title('Template','Fontsize',14);

%% Feature-based alignment
% Extract SURF features and match them
[d1, l1]=iat_surf(img,128);
[d2, l2]=iat_surf(tmp,128);

[map, matches, imgInd, tmpInd]=iat_match_features(d1,d2,.9);
% [map, matches, imgInd, tmpInd]=iat_match_features_mex(d1,d2,.9); %mex-implementation

ptsA=l1(imgInd,1:2);
ptsB=l2(tmpInd,1:2);

% Plot initial correspondences
iat_plot_correspondences(img,tmp,ptsA',ptsB');
title('Initial correspondences','Fontsize',14);

% Apply RANSAC to putative correspondences to get the transform
[inliers, ransacWarp]=iat_ransac(iat_homogeneous_coords(ptsB'),iat_homogeneous_coords(ptsA'),...
    transform, 'tol',.05, 'maxInvalidCount', 10);

%% Plot results
% plot the new correspondences
iat_plot_correspondences(img,tmp,ptsA(inliers,1:2)',ptsB(inliers,1:2)');
title('RANSAC-filtered correspondences','Fontsize',14);

% Create mosaic
ransacMosaic = iat_mosaic(tmp,img,ransacWarp);
figure;imshow(uint8(ransacMosaic));title('Mosaic after feature-based alignment (RANSAC-homography)','Fontsize',14);

% Compute the warped image and visualize the error
[wimage, support] = iat_inverse_warping(img, ransacWarp, transform, 1:size(tmp,2),1:size(tmp,1));

% visualize the error
[imdiff, grayerror] = iat_error2gray(tmp,wimage,support);
figure;imshow(grayerror); title('Error of feature-based alignment (RANSAC)', 'Fontsize', 14);


%% Area-based alignment for refinement

clear par;
% Call ECC initialized by Ransac result
par.transform = transform;
par.iterations = 15;
par.initwarp = ransacWarp;%initialization

% run ECC
[ECCWarp]=iat_ecc(img,tmp,par);

% Create new mosaic
ECCMosaic = iat_mosaic(tmp,img,ECCWarp);
figure;imshow(uint8(ECCMosaic));title('Mosaic after RANSAC-homography+ECC','Fontsize',14);

% Compute the warped image and visualize the error
[wimageECC, supportECC] = iat_inverse_warping(img, ECCWarp, transform, 1:size(tmp,2),1:size(tmp,1));

% visualize the error
[imdiffECC, grayerrorECC] = iat_error2gray(tmp,wimageECC,supportECC);
figure;imshow(grayerrorECC); title('Error of area-based alignment (RANSAC+ECC)', 'Fontsize', 14);

