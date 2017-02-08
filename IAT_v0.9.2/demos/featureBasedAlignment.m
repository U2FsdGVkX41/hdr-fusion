clc;
close all;
clear all;

%% Load and show images
% Read initial image <-> template pair
%img=imread('textImage.png');

%tmp=imread('textTemplate.png');
img=imread('G:\HDRImageSynthesis\PatchBasedHDR_MATLAB_v1.1\Scenes\lady\Img1.tif');
tmp=imread('G:\HDRImageSynthesis\PatchBasedHDR_MATLAB_v1.1\Scenes\lady\Img2.tif');

transform = 'affine'; % it can be approximated by euclidean as well

% Plot both of these images
figure;imshow(img);title('Image','Fontsize',14);
figure;imshow(tmp);title('Template','Fontsize',14);

%% Feature-based alignment
% Extract SURF descriptors
[d1, l1]=iat_surf(img);
[d2, l2]=iat_surf(tmp);

% Match keypoints
%[map, matches, imgInd, tmpInd]=iat_match_features(d1,d2,.7);
[map, matches, imgInd, tmpInd]=iat_match_features_mex(d1,d2,.7);

ptsA=l1(imgInd,1:2);
ptsB=l2(tmpInd,1:2);

%Least-Square
LSWarp = iat_get_affine(iat_homogeneous_coords(ptsB'),iat_homogeneous_coords(ptsA'));

% RANSAC
[inliers, ransacWarp]=iat_ransac(iat_homogeneous_coords(ptsB'),iat_homogeneous_coords(ptsA'),...
transform,'tol',.05, 'maxInvalidCount', 10);

%% Plot results
% Plot initial correspondences
iat_plot_correspondences(img,tmp,ptsA',ptsB');
title('Initial correspondences','Fontsize',14);

%Plot filtered correspondences
iat_plot_correspondences(img,tmp,ptsA(inliers,:)',ptsB(inliers,:)');
title('RANSAC-filtered correspondences','Fontsize',14);

% Compute the warped image and visualize the error
[wimage, support] = iat_inverse_warping(img, ransacWarp, transform, 1:size(tmp,2),1:size(tmp,1));

% Compute the warped image and visualize the error
[wimageLS, supportLS] = iat_inverse_warping(img, LSWarp, transform, 1:size(tmp,2),1:size(tmp,1));

% 
% 
% %plot the warped image
 figure;imshow(uint8(wimage)); title('Warped image by feature-based alignment (RANSAC)', 'Fontsize', 14);
% 
% figure;imshow(uint8(wimageLS)); title('Warped image by feature-based alignment (Least-Squares)', 'Fontsize', 14);
% 
% 
% % visualize the error
% [imdiff, grayerror] = iat_error2gray(tmp,wimage,support);
% figure;imshow(grayerror); title('Error of RANSAC feature-based alignment ', 'Fontsize', 14);
% 
% % visualize the error
% [~, grayerrorLS] = iat_error2gray(tmp,wimage,supportLS);
% figure;imshow(grayerrorLS); title('Error of Least-squares feature-based alignment ', 'Fontsize', 14);
% 
% 
% 
% 
% %% ECC without initialization
% clear par;
% par.levels = 3;
% par.iterations = 25;
% par.transform = transform;
% 
% [ECCWarp]=iat_ecc(img,tmp,par);
% 
% % Compute the warped image and visualize the error
% [wimageECC, supportECC] = iat_inverse_warping(img, ECCWarp, par.transform, 1:size(tmp,2),1:size(tmp,1));
% 
% 
% %plot the warped image
% figure;imshow(uint8(wimageECC)); title('Warped image by ECC', 'Fontsize', 14);
% 
% % visualize the error
% [~, grayerrorECC] = iat_error2gray(tmp,wimageECC,supportECC);
% figure;imshow(grayerrorECC); title('Error of ECC alignment', 'Fontsize', 14);
% 
% 
% 
% %% ECC with rough initialization
% par.levels = 1;
% par.iterations = 25;
% par.transform = transform;
% theta = -pi/180*8;
% par.initwarp = [cos(theta) -sin(theta) -25; sin(theta) cos(theta) 20];
% 
% 
% % warped-image and error due to initialization
% [wimageInit, supportInit] = iat_inverse_warping(img, par.initwarp, par.transform, 1:size(tmp,2),1:size(tmp,1));
% 
% %plot the warped image
% figure;imshow(uint8(wimageInit)); title('Warped image after manual initialization', 'Fontsize', 14);
% 
% % visualize the error
% [~, grayerrorInit] = iat_error2gray(tmp,wimageInit,supportInit);
% figure;imshow(grayerrorInit); title('Error after initialization', 'Fontsize', 14);
% 
% 
% % Run ECC
% [ECCInitWarp]=iat_ecc(img,tmp,par);
% 
% % Compute the warped image and visualize the error
% [wimageECCInit, supportECCInit] = iat_inverse_warping(img, ECCInitWarp, par.transform, 1:size(tmp,2),1:size(tmp,1));
% 
% %plot the warped image
% figure;imshow(uint8(wimageECCInit)); title('Warped image after Initialization+ECC', 'Fontsize', 14);
% 
% % visualize the error
% [~, grayerrorECCInit] = iat_error2gray(tmp,wimageECCInit,supportECCInit);
% figure;imshow(grayerrorECCInit); title('Error of Initialization+ECC alignment', 'Fontsize', 14);
% 
% 
