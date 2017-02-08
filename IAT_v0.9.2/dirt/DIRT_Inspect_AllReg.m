function DIRT_Inspect_AllReg(allReg, varargin)
% DIRT_Inspect_AllReg(allReg, varargin)
%
% Do something with the structure containing the computed registrations
% through the iterations.
%
% Mandatory inputs:
%  - allReg [structure]
%   The structure containing the registrations, as returned by
%   DIRT_Registration in reg.allReg
%
% Optional inputs (properties):
%  - BackupD [string]
%   Backup the difference images using the format string -- an fprintf with
%   the iteration number is used, and the image extension, e.g. png must be
%   included
%  - BackupW [string]
%   Backup the warped images using the format string -- an fprintf with
%   the iteration number is used, and the image extension, e.g. png must be
%   included
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

z = 1;
while z <= length(varargin)
    switch(varargin{z})
        case 'BackupD'
            fmtfn = varargin{z+1};
            for t = 1:length(allReg)
                imwrite(uint8(allReg{t}.D), sprintf(fmtfn, t));
            end;
            z = z + 2;
        case 'BackupW'
            fmtfn = varargin{z+1};
            for t = 1:length(allReg)
                imwrite(uint8(allReg{t}.W), sprintf(fmtfn, t));
            end;
            z = z + 2;
        otherwise
            error(['unknown property' varargin{z}]);
    end;    
end;
