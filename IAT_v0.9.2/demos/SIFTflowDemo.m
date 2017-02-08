% load the imahes
%im1=imread('Mars1HalfSize.png');
%im2=imread('Mars2HalfSize.png');

im1=imread('G:\HDRImageSynthesis\PatchBasedHDR_MATLAB_v1.1\Scenes\lady\Img1.tif');
im2=imread('G:\HDRImageSynthesis\PatchBasedHDR_MATLAB_v1.1\Scenes\lady\Img2.tif');

% Compute SIFT images
patchsize=8; %half of the window size for computing SIFT
gridspacing=1; %sampling step (use 1 for registration)


% compute dense SIFT images
Sift1=iat_dense_sift(im2double(im1),patchsize,gridspacing);
Sift2=iat_dense_sift(im2double(im2),patchsize,gridspacing);

% visualize SIFT image
figure;imshow(iat_sift2rgb(Sift1));title('SIFT image 1');
figure;imshow(iat_sift2rgb(Sift2));title('SIFT image 2');


% SIFT-flow parameters
SIFTflowpara.alpha=2;
SIFTflowpara.d=40;
SIFTflowpara.gamma=0.005;
SIFTflowpara.nlevels=4;
SIFTflowpara.wsize=5;
SIFTflowpara.topwsize=20;
SIFTflowpara.nIterations=60;

% Run the algorithm
tic;[vx,vy,energylist]=iat_SIFTflow(Sift1,Sift2,SIFTflowpara);toc


% VISUALIZE RESULTS

% Keep the pixels that are present in SIFT images 
if gridspacing==1
    Im1=im1(patchsize/2+1:end-patchsize/2,patchsize/2+1:end-patchsize/2,:);
    Im2=im2(patchsize/2+1:end-patchsize/2,patchsize/2+1:end-patchsize/2,:);
else
  im1filt=imfilter(im1,fspecial('gaussian',7,1.),'same','replicate');
  Im1 = im1filt(patchsize/2:gridspacing:end-patchsize/2,patchsize/2:gridspacing:end-patchsize/2,:);
  im2filt=imfilter(im2,fspecial('gaussian',7,1.),'same','replicate');
  Im2 = im2filt(patchsize/2:gridspacing:end-patchsize/2,patchsize/2:gridspacing:end-patchsize/2,:);
end

% warp the image (inverse warping of Im2)
[warpI2, support] = iat_pixel_warping(Im2,vx,vy);

figure;imshow(Im1);title('Image 1');
figure;imshow(uint8(warpI2));title('Warped Image 2');

% visualize alignment error
[~, grayerror] = iat_error2gray(Im1,warpI2,support);
figure;imshow(grayerror);title('Registration error');

% display flow
figure;imshow(iat_flow2rgb(vx,vy));title('SIFT flow field');

