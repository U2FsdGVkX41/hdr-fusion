clc;
clear all;
close all;

%% load ans show images
% Read initial image <-> template pair
img=imread('triangleImage.png');
tmp=imread('triangleTemplate.png');

% Plot both of these images
figure;imshow(img);title('Image', 'Fontsize', 14);
figure;imshow(tmp);title('Template', 'Fontsize', 14);

% choose one between euclidean and affine (the real transform is euclidean)
transform = 'affine';
%transform = 'euclidean';

%% Area-based alignment with ECC
par = [];
par.levels = 1;
par.iterations = 20;
par.transform = transform;

[ECCWarp]=iat_ecc(img,tmp,par);

% Compute the warped image and visualize the error
[wimageECC, supportECC] = iat_inverse_warping(img, ECCWarp, par.transform, 1:size(tmp,2),1:size(tmp,1));


%plot the warped image
figure;imshow(uint8(wimageECC)); title('Warped image by ECC', 'Fontsize', 14);

% visualize the error
[imdiffECC, grayerrorECC] = iat_error2gray(tmp,wimageECC,supportECC);
figure;imshow(grayerrorECC); title('Error of ECC alignment', 'Fontsize', 14);


%% Area-based alignment with Lucas-Kanade
% same parameters
[LKWarp]=iat_LucasKanade(img,tmp,par);

% Compute the warped image and visualize the error
[wimageLK, supportLK] = iat_inverse_warping(img, LKWarp, par.transform, 1:size(tmp,2),1:size(tmp,1));

%plot the warped image
figure;imshow(uint8(wimageLK)); title('Warped image by Lucas-Kanade', 'Fontsize', 14);

% visualize the error
[imdiffLK, grayerrorLK] = iat_error2gray(tmp,wimageLK,supportLK);
figure;imshow(grayerrorLK); title('Error of Lucas-Kanade alignment', 'Fontsize', 14);




