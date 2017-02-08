
addpath(genpath('./GraphSeg/GraphSeg'));
addpath(genpath('./gco-v3.0'));


global resultPath;
global imagesPath;
global imagesExtention;
global handheld;
global refIndex;

%% settings
imagesExtention = 'tif';
evExtention = 'txt';
refIndex = 3;
handheld = 0; % 0 for static camera dataset

if ~isdir(resultPath)
    mkdir(resultPath);
end

imagesPath = 'G:\HDRImageSynthesis\imgdata\PuppetsSequence\';

%% Read images
imagesFileName = dir([imagesPath '*.' imagesExtention]);
imageInfo = imfinfo([imagesPath imagesFileName(1).name]);
imageData = uint8(zeros(imageInfo.Height, imageInfo.Width, 3, size(imagesFileName, 1)));
for imageIndex = 1:size(imagesFileName, 1)
    imageData(:, :, :, imageIndex) = imread([imagesPath imagesFileName(imageIndex).name]);
end

%% Read exposure time
evFileName = dir([imagesPath '*.' evExtention]);
fid = fopen(char([imagesPath evFileName(1).name]));
exposureTime = fscanf(fid, '%f', inf);
exposureTime = 1./exposureTime';
fclose(fid);

%% ECC Alignment
if (handheld == 1)
 mkdir([imagesPath,'align\']);   
transform = 'euclidean';

par = [];
par.levels =    2;
par.iterations = 30;
par.transform = transform;

par1 = [];
par1.levels =    2;
par1.iterations = 10;
par1.transform = transform;

for imageIndex = 1:size(imagesFileName, 1)
    if (imageIndex == refIndex)
        imwrite(imageData(:, :, :, imageIndex),[imagesPath,'align\',num2str(imageIndex),'.png']);
        continue;
    end
[ECCWarp]=iat_ecc(imageData(:, :, :, imageIndex),imageData(:, :, :, refIndex),par);
[wimageECC, supportECC] = iat_inverse_warping(imageData(:, :, :, imageIndex), ECCWarp,...
    par1.transform, 1:size(imageData(:, :, :, refIndex),2),1:size(imageData(:, :, :, refIndex),1));
imageData(:, :, :, imageIndex)=uint8(wimageECC);
imwrite(imageData(:, :, :, imageIndex),[imagesPath,'align\',num2str(imageIndex),'.png']);
end;
end;


%% Salience map
if(handheld==0)
    dos(['Saliency.exe ','"',strrep(imagesPath,'\','/'),'" ','"',imagesExtention,'"']);
else
    dos(['Saliency.exe ','"',strrep(imagesPath,'\','/'),'align/" ','"','png','"']);
end

%% Estimate CRF
CRF = EstimateCRF(imageData, exposureTime);

%% Fusion
resultImage = GhostFreeImageSynthesis(imageData, CRF, exposureTime);
imwrite(resultImage, [resultPath 'result.png']);
