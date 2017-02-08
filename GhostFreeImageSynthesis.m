function result = GhostFreeImageSynthesis(imageData, CRF, exposureTime)

global imagesPath;
global refIndex;
global handheld;

%refIndex = ceil(double(numImages) * 0.5);

imageHeight = size(imageData, 1);
imageWidth = size(imageData, 2);
numChannels = size(imageData, 3);
numImages = size(imageData, 4);


%% Salience Import
[H,W,~,~]=size(imageData);

%sal1=zeros(r,c,n);
if(handheld==0)
    salFileName = dir([imagesPath,'sal\' '*.' 'png']);
    for imageIndex = 1:size(salFileName, 1)
        sal(:,:,imageIndex) = im2double(imread([imagesPath,'sal\',salFileName(imageIndex).name]));
    end
else
      salFileName = dir([imagesPath,'align\sal\' '*.' 'png']);
    for imageIndex = 1:size(salFileName, 1)
        sal(:,:,imageIndex) = im2double(imread([imagesPath,'align\sal\',salFileName(imageIndex).name]));
    end
end

% for imageIndex = 1:numImages
%  %   sal1(:,:,imageIndex) = im2double(imread(['C:\Users\ThinkPad\Desktop\DataSet1\Saliency\A\000',num2str(imageIndex),'_RC.png']));
%    sal(:,:,imageIndex) = im2double(imread([imagesPath,'sal\',imagesFileName,'_RC.png']));
% end

%% Detect moving objects
%objMask = MyObjectDetection(imageData, exposureTime, refIndex, CRF, sal);
ostart = tic;
objMask1 = MyObjectDetectionsal(imageData, exposureTime, CRF,sal);
objMask2 = MyObjectDetectionInversesal(imageData, exposureTime,  CRF,sal);
objMask = objMask1|objMask2;

if(handheld==1)
    
    objMask(1:5,:,:,:) = 0;
    objMask(:,1:5,:,:) = 0;
    objMask(imageHeight-5:imageHeight,:,:,:) = 0;
    objMask(:,imageWidth-5:imageWidth,:,:) = 0;  
end

objMask(:,:,1,1) = bwareaopen(objMask1(:,:,1,1),200);
objMask(:,:,2,1) = bwareaopen(objMask1(:,:,2,1),200);
objMask(:,:,3,1) = bwareaopen(objMask1(:,:,3,1),200);

objMask(:,:,1,2) = bwareaopen(objMask(:,:,1,2),200);
objMask(:,:,2,2) = bwareaopen(objMask(:,:,2,2),200);
objMask(:,:,3,2) = bwareaopen(objMask(:,:,3,2),200);

otime = toc(ostart);

%% Mask extension
estart = tic;

refImage = imageData(:, :, :, refIndex);
refImage_I = rgb2gray(im2double(refImage));
refImageFilted_I = CoherenceFilter(refImage_I, struct('T', 5, 'rho', 3,'Scheme','I', 'sigma', 1));
refSeg = graphSeg(refImageFilted_I, 1.0, 200, 2, 0);

objMaskEx = false(size(objMask));
objMaskRaw = objMaskEx;
for imageIndex = 1:numImages
    if imageIndex ~= refIndex
        objMaskTmp = uint8(squeeze(objMask(:, :, :, imageIndex)));
        objMaskSqueeze = sum(objMaskTmp, 3) >= 1;

        refSegMasked = refSeg(objMaskSqueeze);
        refSegMaskedLabels = unique(refSegMasked(:));
        objMaskEx(:, :, imageIndex) = ismember(refSeg, refSegMaskedLabels);
        objMaskRaw(:,:,imageIndex) = objMaskEx(:,:,imageIndex);
        hardB = ~objMaskEx(:, :, imageIndex);
        hardB = imerode(hardB, strel('disk', 5));              
        objMaskEx(:, :, imageIndex) = MaskExtension(refImageFilted_I, objMaskSqueeze, hardB);
        objMaskEx(:, :, imageIndex) = imdilate(objMaskEx(:, :, imageIndex), strel('disk',3));
    %   objMask(:,:,imageIndex)=~hardB;
    end
end
etime=toc(estart);
% %%%% HDR SYNTHESIS %%%%%%%%%%%%%
% weight = zeros(256, 1);
% for h = 1:256
%     if (h-1) > 128
%         weight(h) = 255 - (h - 1);
%     else
%         weight(h) = (h - 1) - 0;
%     end
% end
% 
% load('G:\HDRImageSynthesis\test161013\logE Arch.mat'); 
% 
% for r = 1:numImages
%     if r ~= refIndex
% %         g_map(:,:,r) = Omega2(:,:,r).* g_map(:,:,r) .* g_map(:,:,idx_ref);
%         objMaskEx(:,:,r) = objMaskEx(:,:,r) .* objMaskEx(:,:,refIndex);
%         se = strel('disk',1);
%         objMaskEx(:,:,r) = imdilate(objMaskEx(:,:,r),se);
%         se = strel('disk',2);
%         objMaskEx(:,:,r) = imerode(objMaskEx(:,:,r),se);
%     end
% end
% objMaskEx(:,:,refIndex) = ones(H,W);
% 
% hdr = zeros(H,W,3);
% for idx_color = 1:3
%     numerator = zeros(H,W);
%     denominator = zeros(H,W);
%     for r = 1:numImages
%         w = objMaskEx(:,:,r).*weight(objMaskEx(:,:,idx_color,r) + 1);
%         denominator = denominator + w;
%         numerator = numerator + w .* exp(logE(:,:,idx_color,r));
%     end
%     Ei = ones(H,W);
%     
%     idx = find(denominator>0);
%     Ei(idx) = numerator(idx) ./ denominator(idx);
%     
%     hdr(:,:,idx_color) = Ei;
% end
% rgb = tonemap(hdr, 'AdjustLightness', [0.01 0.99], 'AdjustSaturation', 1.75);
% rgb = uint8(((double(rgb)/255).^(1/0.75))*255);
%     
%  imwrite(rgb, sprintf('hdr_Arch_rm.png'));
%  save_exr(hdr, sprintf('hdr_Arch_rm.exr'));
% result = rgb;

% objMaskExgaus(:,:,1)=imfilter(objMaskEx(:,:,1),fspecial('gaussian', [5 5], 2),'same');
% objMaskExgaus(:,:,2)=imfilter(objMaskEx(:,:,2),fspecial('gaussian', [5 5], 2),'same');
% objMaskExgaus(:,:,3)=imfilter(objMaskEx(:,:,3),fspecial('gaussian', [5 5], 2),'same');
% objMaskExgaus(:,:,5)=imfilter(objMaskEx(:,:,5),fspecial('gaussian', [5 5], 2),'same');
% objMaskExgaus(:,:,6)=imfilter(objMaskEx(:,:,6),fspecial('gaussian', [5 5], 2),'same');
% objMaskExgaus(:,:,7)=imfilter(objMaskEx(:,:,7),fspecial('gaussian', [5 5], 2),'same');

%% Fusion: part1
fstart = tic;
objMaskEx2(:,:,1)=objMaskEx(:,:,1);
objMaskEx2(:,:,2)=objMaskEx(:,:,1);
for imageIndex = 2:numImages-1
    objMaskEx2(:,:,3*imageIndex-3)=objMaskEx(:,:,imageIndex);
    objMaskEx2(:,:,3*imageIndex-2)=objMaskEx(:,:,imageIndex);
    objMaskEx2(:,:,3*imageIndex-1)=objMaskEx(:,:,imageIndex);
end
objMaskEx2(:,:,3*numImages-3) = objMaskEx(:,:,numImages);
objMaskEx2(:,:,3*numImages-2) = objMaskEx(:,:,numImages);

%% Fusion: part2

imageDataEx(:,:,:,1) = imageData(:,:,:,1);
refLogE = zeros(imageHeight, imageWidth, numChannels);
for colorIndex = 1:numChannels
    CRFIndex = imageData(:, :, colorIndex, 1) + 1;
    refLogE(:, :, colorIndex) = reshape(CRF(CRFIndex, colorIndex), imageHeight, imageWidth) ...
                                - log(exposureTime(1));
    imageDataEx(:,:,colorIndex,2)=...
        double( iCRF( CRF(:,colorIndex), exp(refLogE(:, :, colorIndex)), exposureTime(2) ) );                      
end

for imageIndex = 2:numImages-1
    for colorIndex = 1:numChannels
    CRFIndex = imageData(:, :, colorIndex, imageIndex) + 1;
    refLogE(:, :, colorIndex) = reshape(CRF(CRFIndex, colorIndex), imageHeight, imageWidth) ...
                                - log(exposureTime(imageIndex));
    imageDataEx(:,:,colorIndex,3*imageIndex-3)=...
        double( iCRF( CRF(:,colorIndex), exp(refLogE(:, :, colorIndex)), exposureTime(imageIndex-1)));   
    imageDataEx(:,:,colorIndex,3*imageIndex-2)=...
        double( iCRF( CRF(:,colorIndex), exp(refLogE(:, :, colorIndex)), exposureTime(imageIndex)));
    imageDataEx(:,:,colorIndex,3*imageIndex-1)=...
        double( iCRF( CRF(:,colorIndex), exp(refLogE(:, :, colorIndex)), exposureTime(imageIndex+1)));
    end  
end

for colorIndex = 1:numChannels
    CRFIndex = imageData(:, :, colorIndex, numImages) + 1;
    refLogE(:, :, colorIndex) = reshape(CRF(CRFIndex, colorIndex), imageHeight, imageWidth) ...
                                - log(exposureTime(numImages));
    imageDataEx(:,:,colorIndex,3*numImages-3)=...
        double( iCRF( CRF(:,colorIndex), exp(refLogE(:, :, colorIndex)), exposureTime(numImages-1) ) );                      
end
imageDataEx(:,:,:,3*numImages-2) = imageData(:,:,:,numImages);

result = ExposureFusion(imageDataEx, objMaskEx2);
ftime = toc(fstart);
pause(0.1);
end
