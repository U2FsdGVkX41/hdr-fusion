function [Gr, Gc] = DIRT_GModel_GradientVectors(gmodel, R, C, gopt)
% [Gr, Gc] = DIRT_HomographyGradientVectors(gmodel, R, C, gopt)
%
% Computes the gradient of the geometric motion model on the row and on
% the columns for a given set of points with rows R and columns C, around
% the identity.
%
% Inputs:
%  - R [column n-vector]
%       The rows of the points of interest
%  - C [column n-vector]
%       The columns of the points of interest
%
% Outputs:
%  - Gr [n x ng matrix]
%       The partial derivatives with respect to the rows
%       with ng the number of parameters of the geometric model
%  - Gc [n x ng matrix]
%       The partial derivatives with respect to the columns
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

n = length(R);
switch gmodel
    case 'Homography'
        Gr = [ R C ones(n, 1) zeros(n, 3) -R.^2 -R.*C ];
        Gc = [ zeros(n, 3) R C ones(n, 1) -R.*C -C.^2 ];
    case '3PtHomography'
        Gr = [ R-gopt.sq(1) C-gopt.sq(2) zeros(n, 2) -(R-gopt.sq(1)).^2 -(R-gopt.sq(1)).*(C-gopt.sq(2)) ];
        Gc = [ zeros(n, 2) R-gopt.sq(1) C-gopt.sq(2) -(R-gopt.sq(1)).*(C-gopt.sq(2)) -(C-gopt.sq(2)).^2 ];
    case 'Affine'
        Gr = [ R C ones(n, 1) zeros(n, 3) ];
        Gc = [ zeros(n, 3) R C ones(n, 1) ];
    case 'Rt' % the parameters for the local model are [ angle tr tc ]
        Gr = [ C ones(n, 1) zeros(n, 1) ];
        Gc = [ -R zeros(n, 1) ones(n, 1) ];
    otherwise
        error(['unknown geometric model ' gmodel]);
end;

