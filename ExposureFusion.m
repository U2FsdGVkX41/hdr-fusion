function result = ExposureFusion(imageData,  objMask)

h = [0 1 0; 1 -4 1; 0 1 0];
sig = 0.2;
global refIndex;
imageHeight = size(imageData, 1);
imageWidth = size(imageData, 2);
numImages = size(imageData, 4);
weight = zeros(imageHeight, imageWidth, numImages);
weightSum = zeros(imageHeight, imageWidth);

image = im2double(imageData);

%% Compute weight
for imageIndex = 1:numImages
    imageSingle = image(:, :, :, imageIndex);
    imageSingle_I = rgb2gray(imageSingle);

    C = abs(imfilter(imageSingle_I, h, 'replicate'));
    S = std(imageSingle, 1, 3);
    E_r = exp( -0.5*( (imageSingle(:, :, 1)) - 0.5).^2 / sig.^2 );
    E_g = exp( -0.5*( (imageSingle(:, :, 2)) - 0.5).^2 / sig.^2 );
    E_b = exp( -0.5*( (imageSingle(:, :, 3)) - 0.5).^2 / sig.^2 );
 %   E_sal = exp( -0.5*( (imageSingle(:, :, 3)) - 0.5).^2 / sig.^2 );
    E = E_r .* E_g .* E_b;
    
    weight(:, :, imageIndex) = (C.^1) .* (S.^1) .* (E.^1) .* ~objMask(:,:,imageIndex);%.* (10*(salience(:,:,imageIndex)+1))%.*(10*salience(:,:,imageIndex)+1);
    weightSum = weightSum + weight(:, :, imageIndex);
end
weightSum = weightSum + 1e-12;

for imageIndex = 1:numImages
    weight(:, :, imageIndex) = weight(:, :, imageIndex) ./ weightSum;
end

% for imageIndex = 1:numImages
%    weight(:, :, imageIndex) = weight(:, :, imageIndex).*(salience(:,:,imageIndex)+1);
%    FlattenedData = weight(:, :, imageIndex); 
%    MappedFlattened = mapminmax(FlattenedData, 0, 1); 
%    weight(:, :, imageIndex) = reshape(MappedFlattened, size( weight(:, :, imageIndex))); 
% end

zeroCount = zeros(imageHeight, imageWidth, 'single');
for imageIndex = 1:numImages
    zeroCount = zeroCount + (weight(:, :, imageIndex) <= 0.01);
end

weightResetMask = (zeroCount == numImages);
for imageIndex = 1:numImages
    weightTmp = weight(:, :, imageIndex);
  %  if(imageIndex ~= refIndex)
   if(imageIndex ~= 3*refIndex-3 ||imageIndex ~= 3*refIndex-2 ||imageIndex ~= 3*refIndex-1)
        weightTmp( weightResetMask ) = 0;
    else
        weightTmp( weightResetMask ) = 1;
    end
    weight(:, :, imageIndex) = weightTmp;
end

% for imageIndex = 1:numImages
%    weight(:, :, imageIndex) = weight(:, :, imageIndex).*(salience(:,:,imageIndex)+1);
%    FlattenedData = weight(:, :, imageIndex); 
%    MappedFlattened = mapminmax(FlattenedData, 0, 1); 
%    weight(:, :, imageIndex) = reshape(MappedFlattened, size( weight(:, :, imageIndex))); 
% end

%% Pyramid decomposition
numLevels = floor(log(min(imageHeight, imageWidth)) / log(2)) - 1;
weightPyramid = cell(numImages,1);
imagePyramid = cell(numImages,1);
for imageIndex = 1:numImages
    weightPyramid{imageIndex, 1} = gaussian_pyramid(weight(:, :, imageIndex), numLevels);
    imagePyramid{imageIndex, 1} = laplacian_pyramid(image(:, :, :,imageIndex), numLevels);
end

resultPyramid = laplacian_pyramid(zeros(imageHeight, imageWidth, 3), numLevels);
for imageIndex = 1:numImages
    for levelIndex = 1:numLevels
        for channelIndex = 1:3
            resultPyramid{levelIndex, 1}(:, :, channelIndex) = resultPyramid{levelIndex, 1}(:, :, channelIndex)...
                + weightPyramid{imageIndex, 1}{levelIndex, 1} .* imagePyramid{imageIndex, 1}{levelIndex, 1}(:, :, channelIndex);
        end
    end
end
result = reconstruct_laplacian_pyramid(resultPyramid);
end