
% this demo shows the difference in running time between forwards additive
% and inverse compositional schemes for both Lucas-Kanade and ECC
% algorithms. Times are shown for homography, affine and euclidean cases.


clc;
close all;
clear all;

% Read an image pair
img=imread('chessImage.png');
tmp=imread('chessTemplate.png');

% parameters for ECC and Lucas-Kanade 
par = [];
par.levels =    1;
par.iterations = 50;
par.transform = 'homography';

disp('--------------------Homography--------------------');

%% Lucas-Kanade algorithm FA
tic;[LKWarpFA]=iat_LucasKanade(img,tmp,par); timeLKFA = toc;
disp(['Forwards additive Lucas-Kanade: ' num2str(round(1000*timeLKFA*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');
%% ECC algorithm FA
tic;[ECCWarpFA]=iat_ecc(img,tmp,par); timeECCFA = toc;
disp(['Forwards additive ECC: ' num2str(round(1000*timeECCFA*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');
%% Lucas-Kanade algorithm IC
tic;[LKWarpIC]=iat_LucasKanadeIC(img,tmp,par); timeLKIC = toc;
disp(['Inverse compositional Lucas-Kanade: ' num2str(round(1000*timeLKIC*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');
%% ECC algorithm IC
tic;[ECCWarpIC]=iat_eccIC(img,tmp,par); timeECCIC = toc;
disp(['Inverse compositional ECC: ' num2str(round(1000*timeECCIC*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');

homoTimes = [timeECCIC timeLKIC timeECCFA timeLKFA];


% parameters for ECC and Lucas-Kanade 
par = [];
par.levels =    1;
par.iterations = 50;
par.transform = 'affine';

disp('--------------------Affine--------------------');

%% Lucas-Kanade algorithm FA
tic;[LKWarpFA]=iat_LucasKanade(img,tmp,par); timeLKFA = toc;
disp(['Forwards additive Lucas-Kanade: ' num2str(round(1000*timeLKFA*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');
%% ECC algorithm FA
tic;[ECCWarpFA]=iat_ecc(img,tmp,par); timeECCFA = toc;
disp(['Forwards additive ECC: ' num2str(round(1000*timeECCFA*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');
%% Lucas-Kanade algorithm IC
tic;[LKWarpIC]=iat_LucasKanadeIC(img,tmp,par); timeLKIC = toc;
disp(['Inverse compositional Lucas-Kanade: ' num2str(round(1000*timeLKIC*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');
%% ECC algorithm IC
tic;[ECCWarpIC]=iat_eccIC(img,tmp,par); timeECCIC = toc;
disp(['Inverse compositional ECC: ' num2str(round(1000*timeECCIC*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');

affineTimes = [timeECCIC timeLKIC timeECCFA timeLKFA];

% parameters for ECC and Lucas-Kanade 
par = [];
par.levels =    1;
par.iterations = 50;
par.transform = 'euclidean';

disp('--------------------Euclidean--------------------');

%% Lucas-Kanade algorithm FA
tic;[LKWarpFA]=iat_LucasKanade(img,tmp,par); timeLKFA = toc;
disp(['Forwards additive Lucas-Kanade: ' num2str(round(1000*timeLKFA*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');
%% ECC algorithm FA
tic;[ECCWarpFA]=iat_ecc(img,tmp,par); timeECCFA = toc;
disp(['Forwards additive ECC: ' num2str(round(1000*timeECCFA*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');
%% Lucas-Kanade algorithm IC
tic;[LKWarpIC]=iat_LucasKanadeIC(img,tmp,par); timeLKIC = toc;
disp(['Inverse compositional Lucas-Kanade: ' num2str(round(1000*timeLKIC*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');
%% ECC algorithm IC
tic;[ECCWarpIC]=iat_eccIC(img,tmp,par); timeECCIC = toc;
disp(['Inverse compositional ECC: ' num2str(round(1000*timeECCIC*1000/par.iterations)/1000) 'ms per iteration'])
disp(' ');

euclideanTimes = [timeECCIC timeLKIC timeECCFA timeLKFA];

figure;barh([1,2,3],[euclideanTimes; affineTimes; homoTimes], 'grouped');
set(gca, 'YTickLabel', ['  euclidean ';'    affine  ';'  homography']);
set(gca, 'FontSize', 16);
legend('ECC-IC','LK-IC','ECC-FA','LK-FA','Location','SouthEast');
xlabel('seconds');
title('Running times of 50 iterations');

