function iat_mex()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IAT_MEX() compiles C/C++ files specified within the function. 
% The MEX files are stored in <iatRoot>/mex/<ARCH> folder, 
% where ARCH is what the matlab function 'computer' retuns, i.e.
% ARCH = computer('arch'). E.g., when Linux 64bit is used, the directory is
% <iatRoot>/mex/glnxa64/
% Appropriate message verifies the successful compilation. Otheriwse, the 
% execution breaks and there will appear errors from compilation process.
% You should see 'DONE' after each source filename when being compiled.
%
% -------------------
% Authors: Georgios Evangelidis
% Copyright (C) 2013 Georgios Evangelidis
% All rights reserved.
%
% For any bugs, please contact <georgios@iatool.net> 
%
% This file is part of the IAT library and is made available under
% the terms of the GNU license (see the COPYING file).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=mfilename('fullpath');
rootDir=fileparts(p);

mexDir = [rootDir filesep 'mex' filesep computer('arch') filesep];

if ~exist(mexDir,'dir')
    eval(['mkdir ' mexDir]);
end

% iat_match_features
fprintf('Compiling iat_match_features_mex.cpp......\n');
command = ['mex ' rootDir filesep 'mex' filesep 'iat_match_features_mex.cpp -outdir ' mexDir];
eval(command);
fprintf('DONE\n');

% SIFTflow
fprintf('Compiling SIFTflow.cpp......\n');
SFp = [rootDir filesep 'nonrigid' filesep 'siftflow' filesep 'cpp' filesep]; 
command = ['mex ' SFp 'SIFTflow.cpp ' SFp 'BPFlow.cpp ' SFp 'Stochastic.cpp -outdir ' mexDir ' -output iat_SIFTflow_mex'];
eval(command);
fprintf('DONE\n');

% SimpleFlow 
% to be filled


disp(['MEX functions have been successfully compiled. MEX files have been stored in ' mexDir]);