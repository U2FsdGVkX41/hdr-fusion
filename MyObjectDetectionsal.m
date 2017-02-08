function [ mask ] = MyObjectDetection2( imageData, exposureT, CRF,salience, optParams)

global refIndex;
%% Set default params
if ~exist('optParams', 'var')
    optParams = [];
end
if ~isfield(optParams, 'amplify')
    optParams.amplify = 10.0;
end
if ~isfield(optParams, 'gamma')
    optParams.gamma = 40.0;
end
if ~isfield(optParams, 'threshold')
    optParams.threshold = 5.0;
end
if ~isfield(optParams, 'beta')
    optParams.beta = 1.9;
end

imageHeight = size(imageData, 1);
imageWidth = size(imageData, 2);
numChannels = size(imageData, 3);
numImages = size(imageData, 4);


%% Adjust reference


refLogE = zeros(imageHeight, imageWidth, numChannels);
for colorIndex = 1:numChannels
    CRFIndex = imageData(:, :, colorIndex, refIndex) + 1;
    refLogE(:, :, colorIndex) = reshape(CRF(CRFIndex, colorIndex), imageHeight, imageWidth) ...
                                - log(exposureT(refIndex));
end

refAdjust = zeros(imageHeight, imageWidth, numChannels, numImages);
imageData_g = zeros(imageHeight, imageWidth, numChannels, numImages);
for imageIndex = 1:numImages
    for colorIndex = 1:numChannels
        refAdjust(:, :, colorIndex, imageIndex)...
            = double( iCRF( CRF(:,colorIndex), exp(refLogE(:, :, colorIndex)), exposureT(imageIndex) ) );
    end
      refAdjust(:,:,:,imageIndex) = imfilter(refAdjust(:,:,:,imageIndex),fspecial('gaussian', [3 3],5),'same');
      imageData_g(:,:,:,imageIndex) = imfilter(imageData(:,:,:,imageIndex),fspecial('gaussian', [3 3],5),'same');
   
end

refImage = imageData_g(:, :, :, refIndex);
refArea_l = refImage <= 1;%1;

refArea_h = refImage >= 254;%254;
refArea_w = ~(refArea_l | refArea_h);

%% Dist
dist_l = zeros(size(imageData_g));
dist_w = zeros(size(imageData_g));
dist_h = zeros(size(imageData_g));

for imageIndex = 1:numImages
    if imageIndex == refIndex
        continue;
    end
    
  
    
    dist_w(:, :, :, imageIndex)...
        = abs(refAdjust(:, :, :, imageIndex) - double(imageData_g(:, :, :, imageIndex)));
    
    if imageIndex < refIndex
        dist_l(:, :, :, imageIndex)...
            = abs(imageData_g(:, :, :, imageIndex));
        
        dist_h(:, :, :, imageIndex)...
            = refAdjust(:, :, :, imageIndex) - double(imageData_g(:, :, :, imageIndex));
        dist_h(:, :, :, imageIndex)...
            = max(dist_h(:, :, :, imageIndex), zeros(size(dist_h(:, :, :, imageIndex))));
    else
        dist_l(:, :, :, imageIndex)...
            = double(imageData_g(:, :, :, imageIndex)) - refAdjust(:, :, :, imageIndex);
        dist_l(:, :, :, imageIndex)...
            = max(dist_l(:, :, :, imageIndex), zeros(size(dist_l(:, :, :, imageIndex))));
        
        dist_h(:, :, :, imageIndex)...
            =  abs(255 - imageData_g(:, :, :, imageIndex));
    end
end

dist = zeros(size(dist_l));
refArea_lex = repmat(refArea_l, 1, 1, 1, numImages);
refArea_wex = repmat(refArea_w, 1, 1, 1, numImages);
refArea_hex = repmat(refArea_h, 1, 1, 1, numImages);

dist(refArea_lex) = dist_l(refArea_lex);
dist(refArea_wex) = dist_w(refArea_wex);
dist(refArea_hex) = dist_h(refArea_hex);
% for imageIndex = 1:numImages
%        if imageIndex == refIndex
%         continue;
%        end
%        for dim = 1:3
%         dist_l(:,:,dim,imageIndex) = dist_l(:, :, dim,imageIndex)-10*abs(salience(:,:,imageIndex)-salience(:,:,3));
%         dist_h(:,:,dim,imageIndex) = dist_h(:, :, dim,imageIndex)-10*abs(salience(:,:,imageIndex)-salience(:,:,3));
%         dist_w(:,:,dim,imageIndex) = dist_w(:, :, dim,imageIndex)+abs(salience(:,:,imageIndex)-salience(:,:,3));
%        end
% end

% for imageIndex = 1:numImages
%        if imageIndex == refIndex
%         continue;
%        end
%        [x,y,~]=size(salience(:,:,imageIndex));
%        rows = reshape(salience(:,:,imageIndex),x*y,1);
%        img_mean = mean(rows);
%        salience(:,:,imageIndex)=salience(:,:,imageIndex)/img_mean;
%        for dim = 1:3
%          dist_l(:,:,dim,imageIndex) = dist_l(:, :, dim,imageIndex).*salience(:,:,imageIndex);
% %           FlattenedData = dist(:,:,dim,imageIndex) 
% %           MappedFlattened = mapminmax(FlattenedData, 0, 255); 
% %           dist(:,:,dim,imageIndex) = reshape(MappedFlattened, size( dist(:,:,dim,imageIndex))); 
%        end
% end

%     sal_l = imageData_g(:,:,:,imageIndex)<=25;
%     sal_h = imageData_g(:,:,:,imageIndex)>=230;
%     sal_w = ~(sal_l | sal_h); 

% for imageIndex = 1:numImages
%        if imageIndex == refIndex
%         continue;
%        end
% %        [x,y,~]=size(salience(:,:,imageIndex));
% %        rows = reshape(salience(:,:,imageIndex),x*y,1);
% %        img_mean = mean(rows);
% %        salience(:,:,imageIndex)=salience(:,:,imageIndex)/img_mean;
%        for dim = 1:3      
%          dist(:, :, dim,imageIndex) = dist(:, :, dim,imageIndex)+10*salience(:,:,imageIndex).*sal_w(:,:,dim);
%          dist(:, :, dim,imageIndex) = dist(:, :, dim,imageIndex)-10*salience(:,:,imageIndex).*sal_l(:,:,dim);
%  %        dist(:, :, dim,imageIndex) = dist(:, :, dim,imageIndex)-10*salience(:,:,imageIndex).*sal_h(:,:,dim);
%        end
% end

% for imageIndex = 1:numImages
%        if imageIndex == refIndex
%         continue;
%        end
%        for dim = 1:3
%          dist_h(:,:,dim,imageIndex) = dist_h(:, :, dim,imageIndex)-salience(:,:,imageIndex)*10;
% %           FlattenedData = dist(:,:,dim,imageIndex) 
% %           MappedFlattened = mapminmax(FlattenedData, 0, 255); 
% %           dist(:,:,dim,imageIndex) = reshape(MappedFlattened, size( dist(:,:,dim,imageIndex))); 
%        end
% end

  sal_l = imageData(:,:,:,imageIndex)<=25;
    sal_h = imageData(:,:,:,imageIndex)>=230;
    sal_w = ~(sal_l | sal_h); 
for imageIndex = 1:numImages
       if imageIndex == refIndex
        continue;
       end

       for dim = 1:3      
         dist(:, :, dim,imageIndex) = dist(:, :, dim,imageIndex)+10*salience(:,:,imageIndex).*sal_w(:,:,dim);
 
      end
end

%dist(refArea_lex) = dist_l(refArea_lex);
%dist(refArea_wex) = dist_w(refArea_wex);
%dist(refArea_hex) = dist_h(refArea_hex);
%% Smooth cost
numSites = imageHeight * imageWidth;
handle = GCO_Create(numSites, 2);

GCO_SetSmoothCost(handle, [0, 1; 1, 0]);
dm1 = [1:numSites, 1:numSites];
dm2 = [2:numSites + 1, (1 + imageHeight):(numSites + imageHeight)];
neighborCost = ones(1, numSites * 2) * optParams.amplify;
sp = sparse(dm1, dm2, neighborCost);
neighbors = sp(1:numSites, 1:numSites);
GCO_SetNeighbors(handle, neighbors);

mask = false(size(imageData_g));

for imageIndex = 1:numImages
    if imageIndex == refIndex
        continue;
    end
    
    for colorIndex = 1:numChannels
        distSingle = squeeze(dist(:, :, colorIndex, imageIndex));
        dataCost_0 = zeros(size(distSingle));
        dataCost_1 = zeros(size(distSingle));
    
        %% low
        
        if imageIndex < refIndex
            th_l = optParams.threshold;
            
            below_l = refArea_l(:, :, colorIndex) & (distSingle <= th_l);
            above_l = refArea_l(:, :, colorIndex) & ~below_l;
            dataCost_1(below_l) = 2 * optParams.gamma;
            cost = 2 * optParams.gamma + abs(distSingle - th_l);
            dataCost_0(above_l) = cost(above_l);
        else
        
            lowDist = distSingle(refArea_l(:, :, colorIndex));
            lowDistStd = std(lowDist(:));
            th_l = optParams.beta * lowDistStd;
            
            below_l = refArea_l(:, :, colorIndex) & (distSingle <= th_l);
            above_l = refArea_l(:, :, colorIndex) & ~below_l;
            dataCost_1(below_l) = 2 * optParams.gamma;
            cost = abs(distSingle - th_l);
            dataCost_0(above_l) = cost(above_l);
        end
    
        %% well
        wellDist = distSingle(refArea_w(:, :, colorIndex));
        wellDistStd = std(wellDist(:));
     %   salienceStd= std(salience(:,:,numImages));
        threshold_w = optParams.beta * wellDistStd;% *salienceStd;
        
        below_w = refArea_w(:, :, colorIndex) & (distSingle <= threshold_w);
        above_w = refArea_w(:, :, colorIndex) & ~below_w;
        dataCost_1(below_w) = 2 * optParams.gamma;
        cost = abs(distSingle - threshold_w);
        dataCost_0(above_w) = cost(above_w);
    
        %% high
        if imageIndex < refIndex
            highDist = distSingle(refArea_h(:, :, colorIndex));
            highDistStd = std(highDist(:));
            th_h = optParams.beta * highDistStd;
            
            below_h = refArea_h(:, :, colorIndex) & (distSingle <= th_h);
            above_h = refArea_h(:, :, colorIndex) & ~below_h;
            dataCost_1(below_h) = 2 * optParams.gamma;
            cost = abs(distSingle - th_h);
            dataCost_0(above_h) = cost(above_h);
            
        else
            th_h = optParams.threshold;
            
            below_h = refArea_h(:, :, colorIndex) & (distSingle <= th_h);
            above_h = refArea_h(:, :, colorIndex) & ~below_h;
            dataCost_1(below_h) = 2 * optParams.gamma;
            cost = 2 * optParams.gamma + abs(distSingle - th_h);
            dataCost_0(above_h) = cost(above_h);
        end
            
    
        %% expansion
        GCO_SetDataCost(handle, (optParams.amplify / optParams.gamma) * [dataCost_0(:), dataCost_1(:)]');
        GCO_Expansion(handle);
        mask(:, :, colorIndex, imageIndex)...
            = reshape(GCO_GetLabeling(handle) == 2, imageHeight, imageWidth);
%         mask(1,:,colorIndex,imageIndex) = 0;
%         mask(:,imageWidth,colorIndex,imageIndex) = 0;
%         mask(:,1,colorIndex,imageIndex) = 0;
%         mask(imageHeight,:,colorIndex,imageIndex) = 0;  
    end
end
pause(0.1);
end

