clc;
close all;
clear all;

%% HOMOGRAPHY alignment 
% Read initial image <-> template pair
img=imread('target.jpg'); % to be warped
tmp=imread('source.jpg'); 


transform = 'homography';

% Plot images
figure;imshow(img);title('Image', 'Fontsize', 14);
figure;imshow(tmp);title('Template', 'Fontsize', 14);


% parameters
par = [];
par.iterations = 100;
par.transform = transform;
par.photoModel = 'SingleGainAndBias'; 

[warpDIC, RoiDIC] = iat_dic(img, tmp, par);

% Compute the warped image and visualize the error
[wimageDIC, supportDIC] = iat_inverse_warping(img, warpDIC, transform, 1:size(tmp,2),1:size(tmp,1));

%plot the warped image
figure;imshow(uint8(wimageDIC)); title('Warped image after DIC', 'Fontsize', 14);

% visualize the error
[imdiffDIC, grayerrorDIC] = iat_error2gray(tmp,wimageDIC,supportDIC);
figure;imshow(grayerrorDIC); title('Error of DIC alignment (without photometric correction)', 'Fontsize', 14);

[imerr] = iat_error2rgb( uint8(tmp), uint8(wimageDIC));
figure;imshow(imerr); title('Colorized alignment error (without photometric correction)','Fontsize',14);

figure; imshow(RoiDIC); title('ROI: Template pixels used for alignment', 'Fontsize', 14);
 

disp('Press any key to run the second example');
pause;

%% Euclidean alignment (with photometric distortion)

clc;
close all;
clear all;

img=imread('BaboonImage.png');
%tmp=imread('BaboonTemplate.png');
tmp=imread('BaboonTemplateIntensityChange.png');
 
transform = 'euclidean';

% Plot images
figure;imshow(img);title('Image', 'Fontsize', 14);
figure;imshow(tmp);title('Template', 'Fontsize', 14);


% parameters
par = [];
par.iterations = 100;
par.transform = transform;
par.photoModel = 'SingleGainAndBias'; 
par.maskBorder = 20;

[warpDIC, RoiDIC] = iat_dic(img, tmp, par);

% Compute the warped image and visualize the error
[wimageDIC, supportDIC] = iat_inverse_warping(img, warpDIC, transform, 1:size(tmp,2),1:size(tmp,1));

%plot the warped image
figure;imshow(uint8(wimageDIC)); title('Warped image after DIC', 'Fontsize', 14);

% visualize the error
[imdiffDIC, grayerrorDIC] = iat_error2gray(tmp,wimageDIC,supportDIC);
figure;imshow(grayerrorDIC); title('Error of DIC alignment (without photometric correction)', 'Fontsize', 14);

[imerr] = iat_error2rgb( uint8(tmp), uint8(wimageDIC));
figure;imshow(imerr); title('Colorized alignment error (without photometric correction)','Fontsize',14);

figure; imshow(RoiDIC); title('ROI: Template pixels used for alignment', 'Fontsize', 14);
