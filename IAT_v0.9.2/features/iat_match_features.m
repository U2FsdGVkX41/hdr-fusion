function [ mapping, nummatches, IndexA, IndexB ] = iat_match_features( descA, descB, ratio )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [MAPPING, NUMMATCHES, INDEXA, INDEXB] = IAT_MATCH_FEATURES(DESCA, DESCB, RATIO)
% IAT_MATCH_FEATURES implements descriptor matching between two sets of
% descriptors (features). The inverse cosine (arccosine) of the angle between two
% descriptors is used for their matching evaluation. Note that the function
% resolves double matches
%
% -->Input:
% DESCA:                NxM array with N descriptors of length M that
%                       correspond to the first image
% DESCB:                KxM array with K descriptors of length M that
%                       correspond to the second image
% RATIO:                The maximum acceptable ratio between nearest and
%                       second_nearest angle for potential match
%                       (default: 1).
%
% -->Output:
% MAPPING:              A Nx1 vector that defines the mapping between DESCA
%                       and DESCB. The indices of MAPPING correspond to
%                       the indices of DESCA, while the entries of MAPPING
%                       correspond to the indicees of DESCB.
%                       Zero entries correspond to unmatched descriptors.
% NUMMATCHES:           Number of confirmed matches
% INDEXA:               A vector with indices of matched descriptors of DESCA
% INDEXB:               A vector with indices of matched descriptors of DESCB (w.r.t. INDEXA)
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

% 0<=ratio<=1
% Default: ratio=1
if nargin==2
    ratio=1;
elseif nargin==3
    if (ratio<=0) || (ratio>1)
        error('iat_match_descriptors: ratio must be in the interval (0 1]');
    end
end


% Number of descriptors for each image
NoDescA=size(descA,1);
NoDescB=size(descB,1);

% Vector that holds the mapping between the descriptors
mapping=zeros(NoDescA,1);

D = zeros(NoDescA, NoDescB);
HistB = zeros(1,NoDescB);
for i=1:NoDescA
    d=zeros(NoDescB,1);
    for j=1:NoDescB
        % Calculate normalized correlations of each descriptor from the
        % first image with every descriptor from the second image
        d(j)=dot(descA(i,:),descB(j,:))/norm(descB(j,:));
        
    end
    acosd = acos(d/norm(descA(i,:)));
    D(i,:) = acosd';
    [sd1, index1] = min(acosd); %best
    acosd(index1) = inf;
    [sd2, index2] = min(acosd); %second-best
    
    % If best<ratio*second_best => match is accepted
    if (d(index2)>0 && sd1/sd2<ratio)
        HistB(index1) = HistB(index1)+1; %counts of correspondences with index1
        %if HistB(index1)<2
            mapping(i) = index1; % save only single matches here
        %end
    end
end


multiIndB = find(HistB>1); % indices from descB with duplicate matches

% when dublicates do an inverse mapping 
if ~isempty(multiIndB)
    for i = 1:length(multiIndB)
        multiMap = (mapping==multiIndB(i));
        mapping(multiMap)=0;
        [~,newIndexA] = min(D(:,multiIndB(i)));
        if mapping(newIndexA)==0 % check availability
            mapping(newIndexA) = multiIndB(i);

        end
    end
end


IndexA = find(mapping>0);
IndexB = mapping(IndexA);



% Number of successful matches
nummatches=numel(IndexA);

disp(['Feature matching: ' num2str(nummatches) ' correspondences found.'])
if nummatches<10
    disp(['iat_match_descriptors: Increase the ratio to get more correspondences'])
end
end