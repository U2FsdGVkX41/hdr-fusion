function stackOut = GrossbergSampling(dir_name, format, stack, nSamples)
%
%       stackOut = GrossbergSampling(dir_name, format, stack, nSamples)
%
%
%        Input:
%           -dir_name: the folder name where the stack is stored as a
%           series of LDR images.
%           -format: an LDR format for reading LDR images in the current directory 
%           -stack: a stack of LDR images; 4-D array where values are
%           -nSamples: the number of samples for sampling the stack
%
%        Output:
%           -stackOut: a stack of LDR samples for Debevec and Malik method
%           (gsolve.m)
%
%     Copyright (C) 2013  Francesco Banterle
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

if(~exist('stack','var')&&~exist('stack_exposure','var'))
    %Read images from the current directory
    stack = ReadLDRStackHistogram(dir_name, format);
end

if(~exist('nSamples','var'))
    nSamples = 100;
end

if(nSamples<1)
    nSamples = 100;
end

[h_bins, col, stackSize] = size(stack);

%Compute CDF
%figure(4);
%hold on;
for i=1:stackSize
    for j=1:col
        h_cdf = cumsum(stack(:,j,i));
        stack(:,j,i) = h_cdf/max(h_cdf(:));
    end
%    plot(stack(:,1,i));
end
%figure(1);
stackOut = zeros(nSamples, stackSize, col);

u = 0:(1.0/(nSamples-1)):1;

for i=1:nSamples
    for j=1:col
        for k=1:stackSize
           [p, val] = min(abs(stack(:,j,k)-u(i)));
           stackOut(i,k,j) = val-1;
        end
    end
end

end