function [m,n] = size(X)
%size (overloaded)
%
% Author: Alexandre 
[m,n] = size(X.data(1).value);
