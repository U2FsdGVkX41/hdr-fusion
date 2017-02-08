clc;
close all;
clear all;

%% Load and show images
% Read initial image <-> template pair
%img=imread('BaboonImage.png');
%tmp=imread('BaboonTemplate.png');
%tmp=imread('BaboonTemplateIntensityChange.png');

img=imread('G:\161228\tmp\P1000584_DxO.tif');
tmp=imread('G:\161228\tmp\P1000585_DxO.tif');


% Plot both of these images
%figure;imshow(img);title('Image','Fontsize',14);
%figure;imshow(tmp);title('Template','Fontsize',14);

transform = 'euclidean';


% parameters for ECC and Lucas-Kanade 
par = [];
par.levels =    2;
par.iterations = 30;
par.transform = transform;

par1 = [];
par1.levels =    2;
par1.iterations = 10;
par1.transform = transform;

%% Lucas-Kanade algorithm
%[LKWarp]=iat_LucasKanade(img,tmp,par);

% Compute the warped image and visualize the error
%[wimageLK, supportLK] = iat_inverse_warping(img, LKWarp, par.transform, 1:size(tmp,2),1:size(tmp,1));

%plot the warped image
%figure;imshow(uint8(wimageLK)); title('Warped image by Lucas-Kanade', 'Fontsize', 14);

% draw mosaic
%LKMosaic = iat_mosaic(tmp,img,[LKWarp; 0 0 1]);
%figure;imshow(uint8(LKMosaic));title('Mosaic after Lucas-Kanade','Fontsize',14);


%% ECC algorithm
pstart=tic;
[ECCWarp]=iat_ecc(img,tmp,par);

% Compute the warped image and visualize the error
[wimageECC, supportECC] = iat_inverse_warping(img, ECCWarp, par.transform, 1:size(tmp,2),1:size(tmp,1));

% [ECCWarp2]=iat_ecc(img2,tmp,par1);
% [wimageECC2, supportECC2] = iat_inverse_warping(img2, ECCWarp2, par1.transform, 1:size(tmp,2),1:size(tmp,1));
% [ECCWarp3]=iat_ecc(img3,tmp,par1);
% [wimageECC3, supportECC3] = iat_inverse_warping(img3, ECCWarp3, par1.transform, 1:size(tmp,2),1:size(tmp,1));
% [ECCWarp4]=iat_ecc(img4,tmp,par1);
% [wimageECC4, supportECC4] = iat_inverse_warping(img4, ECCWarp4, par1.transform, 1:size(tmp,2),1:size(tmp,1));
% [ECCWarp5]=iat_ecc(img5,tmp,par1);
% [wimageECC5, supportECC5] = iat_inverse_warping(img5, ECCWarp5, par1.transform, 1:size(tmp,2),1:size(tmp,1));
% [ECCWarp6]=iat_ecc(img6,tmp,par1);
% [wimageECC6, supportECC6] = iat_inverse_warping(img6, ECCWarp6, par1.transform, 1:size(tmp,2),1:size(tmp,1));
ptime=toc(pstart);
%plot the warped image

%figure;imshow(uint8(wimageECC)); title('Warped image by ECC', 'Fontsize', 14);

% draw mosaic
%ECCMosaic = iat_mosaic(tmp,img,[ECCWarp; 0 0 1]);
%figure;imshow(uint8(ECCMosaic));title('Mosaic after ECC','Fontsize',14);

