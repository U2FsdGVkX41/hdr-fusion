function stackOut = ComputeLDRStackHistogram(stack)
%
%       stackOut = ComputeLDRStackHistogram(stack)
%
%
%        Input:
%           -dir_name: the folder name where the stack is stored as a
%           series of LDR images.
%           -format: an LDR format for reading LDR images in the current directory 
%
%        Output:
%           -stackOut: a stack of LDR image histograms
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

[r,c,col,n] = size(stack);

stackOut = zeros(256,col,n);

for i=1:n
    %store in the stack
    for j=1:col
        stackOut(:,j,i) = imhist(uint8(stack(:,:,j,i)));
    end
end

end