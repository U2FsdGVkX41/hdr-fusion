function x = vect(X)
% x = vect(X);
%
% Vectorizes (column-wise) the matrix X using x = X(:).
% It is useful for vectorizing eg submatrices, as in x = vect(A(1:3, :)).
%
% This function goes with the DIRT Matlab image registration package.
% Please cite my paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

x = X(:);