function newparam = iat_merge_param(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEWPARAM = IAT_MERGE_PARAM(PARAM1,PARAM2,...)
% IAT_MERGE_PARAM merges parameter sets PARAM1, PARAM2,... into NEWPARAM.
% Each parameter set as well as the new one is a struct with field names. 
% If a field is defined multiple times, the last input argument is 
% taken into account.
%
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

newparam = varargin{1};
if ~isstruct(newparam);
    error('iat_merge_param:1st parameter set is not a structure');
end

for i = 2:length(varargin)
    temp = varargin{i};
    if ~isempty(temp)
        if ~isstruct(temp);
            
            switch i
                case 2
                    error('iat_merge_param:2nd parameter set is not a structure');
                case 3
                    error('iat_merge_param:3rd parameter set is not a structure');
                otherwise
                    error(['iat_merge_param:' num2str(i) '-th parameter set is not a structure']);
            end
        end
        for f = fieldnames(temp)'
            %newparam = setfield(newparam,f{1},getfield(temp,f{1}));
newparam.(f{1}) = temp.(f{1});
        end
    end
end

