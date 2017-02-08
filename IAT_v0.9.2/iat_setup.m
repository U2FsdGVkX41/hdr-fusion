function iat_setup(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IAT_SETUP adds the IAT toolbox to MATLAB path
%
% IAT_SETUP('FOREVER') adds the IAT toolbox permanently
% -------------------
% Authors: Georgios Evangelidis, Panagiotis Anatolitis
% Copyright (C) 2013 Georgios Evangelidis
% All rights reserved.
%
% For any bugs, please contact <georgios.evangelidis@inria.fr> or
% <anatolitis@ceid.upatras.gr>
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=mfilename('fullpath');
rootDir=fileparts(p);

if ~isempty(rootDir)
    addpath(genpath(rootDir));
    disp('Image Alignment Toolbox has been successfully installed!');
end

if nargin>0
    if strcmpi('forever',varargin{1})
     savepath;
    end
end