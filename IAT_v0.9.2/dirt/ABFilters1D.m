function h = ABFilters1D(filter_name, sigma, hsz)
% h = ABFilters1D(filter_name, sigma, hsz)
%
% Generate some common 1D filter masks.
%
% Inputs:
%  - filter_name [string]
%       the filter, amongst:
%          * der_x
%          * der_y
%          * gauss_x
%          * gauss_y
%          * gauss_der_x
%          * gauss_der_y
%  - (opt) sigma [scalar]
%       the variance for Gaussian filters. Default: sigma = 1 pixel
%  - (opt) hsz [scalar]
%       half the size of the filter. Default: hsz = 3*sigma pixels
%
% Outputs:
%  - h [1D array]
%       the filter (row vector for x and column vector for y)
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

if ~exist('sigma', 'var') sigma = 1; end;
if ~exist('hsz', 'var') hsz = max(1, 3*sigma); end;

switch filter_name
    
    case 'der_x'
        h = .5 * [-1 0 1];
        
    case 'der_y'
        h = .5 * [-1 ; 0 ; 1];
        
    case 'gauss_x'
        x = -hsz:hsz;
        h = exp(- (x.^2) / (2*sigma^2));
        h = h / sum(h);
        
    case 'gauss_y'
        x = (-hsz:hsz)';
        h = exp(- (x.^2) / (2*sigma^2));
        h = h / sum(h);

    case 'gauss_der_x'
        x = -hsz:hsz;
        h = x .* exp(- (x.^2) / (2*sigma^2));
        h = h / sum(abs(h));
        
    case 'gauss_der_y'
        x = (-hsz:hsz)';
        h = x .* exp(- (x.^2) / (2*sigma^2));
        h = h / sum(abs(h));        
        
    otherwise
        error(['unknown filter' filter_name]);
    
end;
    