degree = get_degree(X)
%  METHOD: get_degree
%
% Author: Alexandre Felipe
% 2014, Dec, 8
%
% Return the degree of the polynomial.

  degree = cellfun(@sum,X.data(1).exponent)
