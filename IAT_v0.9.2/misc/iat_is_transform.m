function [ check ] = iat_is_transform( input )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CHECK]=iat_is_transform(INPUT)
% IAT_IS_TRANSFORM checks if the given string is a valid transform name
% Valid names: {'translation', 'euclidean', 'similarity', 'affine',
% 'homography'}
%
% -->Input:
% INPUT:                The input string.
%
% -->Output:
% CHECK:                A logical value that shows if the input was a valid
%                       geometric transform string.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~ischar(input))
    error('iat_is_transform: Input is not a string');
end

check = strcmpi(input,'translation')||strcmpi(input,'euclidean') ...
      ||strcmpi(input,'similarity') ||strcmpi(input,'affine') ...
      ||strcmpi(input,'homography');

end

