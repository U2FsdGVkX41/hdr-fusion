% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

% This is an example of how to use the registration package.
% Two images are loaded, the user specifies a region of interest, and the
% images are registered with an homography as geometric model and different 
% photometric models.

% load the images
sI_uint8 = imread('source.jpg');
tI_uint8 = imread('target.jpg');

% convert them to gray-level
sI_uint8_gl = rgb2gray(sI_uint8);
tI_uint8_gl = rgb2gray(tI_uint8);

% let the user enter a polygonal delineation for the ROI
M = DIRT_MaskManual(sI_uint8);

% define the pixels of interest as pixels nearby edges
M = M & DIRT_MaskEdges(sI_uint8_gl, 2);

% register using intensity only - photometric model: identity
pre1 = DIRT_Precomputation(double(sI_uint8_gl), 'ROI', DIRT_Mask2ROI(M), 'no_poptim');
reg1 = DIRT_Registration(double(tI_uint8_gl), pre1);

% register using intensity only - photometric model: gain and bias
pre2 = DIRT_Precomputation(double(sI_uint8_gl), 'ROI', DIRT_Mask2ROI(M), 'pmodel', 'GainAndBias');
reg2 = DIRT_Registration(double(tI_uint8_gl), pre2);

% register using color - photometric model: identity
pre3 = DIRT_Precomputation(double(sI_uint8), 'ROI', DIRT_Mask2ROI(M), 'no_poptim');
reg3 = DIRT_Registration(double(tI_uint8), pre3);

% register using color - photometric model: single gain and bias
pre4 = DIRT_Precomputation(double(sI_uint8), 'ROI', DIRT_Mask2ROI(M), 'pmodel', 'SingleGainAndBias');
reg4 = DIRT_Registration(double(tI_uint8), pre4);

% register using color - photometric model: multiple gains and biases
pre5 = DIRT_Precomputation(double(sI_uint8), 'ROI', DIRT_Mask2ROI(M), 'pmodel', 'MultipleGainsAndBiases');
reg5 = DIRT_Registration(double(tI_uint8), pre5);

% register using color - photometric model: full affine channel mixing
pre6 = DIRT_Precomputation(double(sI_uint8), 'ROI', DIRT_Mask2ROI(M), 'pmodel', 'Affine');
reg6 = DIRT_Registration(double(tI_uint8), pre6);

% print a summary
fprintf('--- Summary ---\n');
fprintf('-- Grey-level methods --\n');
fprintf('  Identity:                   %03d iterations | error at convergence: %06.2f | CPU time: %06.2f seconds\n', reg1.it, reg1.e_roi, reg1.cpu_time);
fprintf('  Gain and bias:              %03d iterations | error at convergence: %06.2f | CPU time: %06.2f seconds\n', reg2.it, reg2.e_roi, reg2.cpu_time);
fprintf('-- Color methods --\n');
fprintf('  Identity:                   %03d iterations | error at convergence: %06.2f | CPU time: %06.2f seconds\n', reg3.it, reg3.e_roi, reg3.cpu_time);
fprintf('  Single gain and bias:       %03d iterations | error at convergence: %06.2f | CPU time: %06.2f seconds\n', reg4.it, reg4.e_roi, reg4.cpu_time);
fprintf('  Multiple gains and biases:  %03d iterations | error at convergence: %06.2f | CPU time: %06.2f seconds\n', reg5.it, reg5.e_roi, reg5.cpu_time);
fprintf('  Full affine channel mixing: %03d iterations | error at convergence: %06.2f | CPU time: %06.2f seconds\n', reg6.it, reg6.e_roi, reg6.cpu_time);
fprintf('---------------\n');
