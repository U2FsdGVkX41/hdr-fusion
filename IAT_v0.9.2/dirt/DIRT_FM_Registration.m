function reg_FM = DIRT_FM_Registration(tI, pre_FM, tR, tC, varargin)
% reg_FM = DIRT_FM_Registration(tI, pre_FM, tR, tC, varargin)
%
% Compute the SSD and OptSSD between the points in the source image in pre_FM 
% and the points (tR,tC) in the target image.
%
% Mandatory inputs:
%  - tI [image]
%       The source image
%  - tR, tC [column vectors]
%       Rows and columns in the target image of the points to be matched
%
% Optional inputs (properties):
%
% Outputs:
%  - reg_FM [struct]
%       Contains the computed RawSSD and OptSSD scores
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

% check number of parameters
if nargin < 4 error('not enough input arguments'); end;

% setup default parameters and parse input parameters
optReg_FM = ParseVarargin(pre_FM.optPre_FM, varargin{:});

m = length(tR);

for sj = 1:pre_FM.m

    if optReg_FM.verb fprintf('sj = %04d / %04d\n', sj, pre_FM.m); end;
    
    for tj = 1:m
        if optReg_FM.verb fprintf('  tj = %04d / %04d\n', tj, m); end;
        
        gopt.tq = [tR(tj);tC(tj)];
        
        reg = DIRT_Registration(tI, pre_FM.pre{sj}, 'verb', 0, 'max_nb_it', 1, 'gopt', gopt);
        reg_FM.RawSSD(sj,tj) = reg.e_roi_init;
        if reg.e_roi_init < optReg_FM.max_RawSSD4OptSSD && norm(pre_FM.pre{sj}.optPre.gopt.sq-gopt.tq) < optReg_FM.max_disparity
            reg = DIRT_Registration(tI, pre_FM.pre{sj}, 'verb', 0, 'min_step_size', optReg_FM.min_step_size, ...
                'max_nb_it', optReg_FM.max_nb_it, 'gopt', gopt);
        end;            
        if reg.err            
            reg_FM.OptSSD(sj,tj) = reg.e_roi_init;
        else
            reg_FM.OptSSD(sj,tj) = reg.e_roi;
        end;            
    end;
    
end;

    
function optReg_FM = ParseVarargin(optPre_FM, varargin)

% default values
optReg_FM.max_RawSSD4OptSSD = 30;
optReg_FM.max_disparity = 300;
optReg_FM.verb = 1;
optReg_FM.min_step_size = 1e-2;
optReg_FM.max_nb_it = 10;

% varargin parsing
z = 1;
while z <= length(varargin)
    switch(lower(varargin{z}))
        case 'max_rawssd4optssd'
            optReg_FM.max_RawSSD4OptSSD = varargin{z+1};
            z = z + 2;
        case 'max_disparity'
            optReg_FM.max_disparity = varargin{z+1};
            z = z + 2;
        case 'verb'
            optReg_FM.verb = varargin{z+1};
            z = z + 2;
        case 'min_step_size'
            optReg_FM.min_step_size = varargin{z+1};
            z = z + 2;
        case 'max_nb_it'
            optReg_FM.max_nb_it = varargin{z+1};
            z = z + 2;
        otherwise
            error(['unknown property ' varargin{z}]);
    end;
end;
