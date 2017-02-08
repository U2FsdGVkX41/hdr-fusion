function [ labeling ] = MaskExtension( image_I, hardMask_o, hardMask_b, optParams )

%% Set default params
if ~exist('optParams', 'var')
    optParams = [];
end
if ~isfield(optParams, 'amplify')
    optParams.amplify = 1e3;
end
if ~isfield(optParams, 'lambda')
    optParams.lambda = 50.0;
end
if ~isfield(optParams, 'gamma')
    optParams.gamma = 10;
end

sigma = std(image_I(:)) * optParams.lambda;

image_I = image_I * 255.0;
imageHeight = size(image_I, 1);
imageWidth = size(image_I, 2);
numSites = numel(image_I);

handle = GCO_Create(numSites, 2);
%% Smooth cost
GCO_SetSmoothCost(handle, [0, 1; 1, 0]);

dm1 = [1:numSites, 1:numSites];
dm2 = [2:numSites + 1, (1 + imageHeight):(numSites + imageHeight)];
values = zeros(1, numSites * 2);

image_d = zeros(size(image_I));
image_d(1:imageHeight-1, :) = image_I(2:imageHeight, :);
values_d = exp(- ( abs(double(image_d) - double(image_I)) .^2.0 ./ (2.0*sigma*sigma) ) );
values(1:numSites) = values_d(:);

image_r = zeros(size(image_I));
image_r(:, 1:imageWidth-1) = image_I(:, 2:imageWidth);
values_r = exp(- ( abs(double(image_r) - double(image_I)) .^2.0 ./ (2.0*sigma*sigma) ) );
values(1 + numSites:numSites * 2) = values_r(:);

sp = sparse(dm1, dm2, round( optParams.amplify * values ));
neighbors = sp(1:numSites, 1:numSites);
GCO_SetNeighbors(handle, neighbors);

%% Data cost
cost_0 = zeros(size(image_I));
cost_0(hardMask_o) = 10000;
cost_0(hardMask_b) = 0;
cost_1 = zeros(size(image_I));
cost_1(hardMask_o) = 0;
cost_1(hardMask_b) = 10000;
GCO_SetDataCost(handle, (optParams.amplify / optParams.gamma) * [cost_0(:), cost_1(:)]');

%% Expansion
GCO_Expansion(handle);

labeling = reshape(GCO_GetLabeling(handle) == 2, size(image_I));

end
