function pre_FM = DIRT_FM_Precomputation(sI, sR, sC, varargin)
% pre_FM = DIRT_FM_Precomputation(sI, sR, sC, varargin)
%
% Do some precomputations for two-image feature matching using the
% DIRT_FM_Registration function.
%
% Mandatory inputs:
%  - sI [image]
%       The source image
%  - sR, sC [column vectors]
%       Rows and columns in the source image of the points to be matched
%
% Optional inputs (properties):
%
% Outputs:
%  - pre_FM [struct]
%       Contains the precomputations
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

% check number of parameters
if nargin < 3, error('not enough input arguments'); end;

% setup default parameters and parse input parameters
pre_FM.optPre_FM = ParseVarargin(varargin{:});

% store the source image data
pre_FM.sI = sI;

% the number of features
pre_FM.m = length(sR);

% % the image derivatives
% if pre.optPre.filter_image_der_sigma == 0
%     filter = ABFilters1D('der_x');
% else
%     filter = ABFilters1D('gauss_der_x', pre.optPre.filter_image_der_sigma);
% end;
% if ~pre.optPre.color
%     pre.sIc = conv2(sI, - filter, 'same');
%     pre.sIr = conv2(sI, - filter', 'same');
% else
%     for c = 1:3
%         pre.sIc(:,:,c) = conv2(sI(:, :, c), - filter, 'same');
%         pre.sIr(:,:,c) = conv2(sI(:, :, c), - filter', 'same');
%     end;    
% end;


for sj = 1:pre_FM.m
    
    if pre_FM.optPre_FM.verb, fprintf('sj = %04d / %04d\n', sj, pre_FM.m); end;
    
    M = DIRT_MaskRectangle(sI, round(sR(sj)), round(sC(sj)), pre_FM.optPre_FM.HalfMaskSize, pre_FM.optPre_FM.HalfMaskSize);
    gopt.sq = [sR(sj);sC(sj)];
lear    pre_FM.pre{sj} = DIRT_Precomputation(sI, DIRT_Mask2ROI(M), 'no_poptim', ...
        'gmodel', pre_FM.optPre_FM.gmodel, ...
        'gopt', gopt, ...
        'verb', 0);
        
end;

function optPre_FM = ParseVarargin(varargin)

% default values
optPre_FM.HalfMaskSize = 5;
optPre_FM.gmodel = '3PtHomography';
optPre_FM.verb = 1;

% varargin parsing
z = 1;
while z <= length(varargin)
    switch(lower(varargin{z}))
        case 'halfmasksize'
            optPre_FM.HalfMaskSize = varargin{z+1};
            z = z + 2;
        case 'gmodel'
            optPre_FM.gmodel = varargin{z+1};
            z = z + 2;
        case 'verb'
            optPre_FM.verb = varargin{z+1};
            z = z + 2;
        otherwise
            error(['unknown property ' varargin{z}]);
    end;
end;
