% ALMOSTEQUAL  Check if two arrays are almost equal
%
%   ALMOSTEQUAL(X, Y) returns true if X and Y are arrays of the same size and
%   with the same elements except for differences that may be attributed to
%   rounding errors. The maximum allowed relative difference is 5 x 10^(-14) and
%   the maximum allowed absolute difference is also 5 x 10^(-14). The relative
%   difference is measured relative to the element in X or Y with largest
%   absolute value.
%
%   ALMOSTEQUAL(X, Y, TOL) uses TOL instead of 5 x 10^(-14)

function close = almostequal(x,y,tol)
  if nargin < 3, tol = 5e-14; end
  if iscell(x)
    if ~iscell(y), close=false; return, end
    x = cell2mat(x);
    y = cell2mat(y);
  end
  M = max([1; abs(x(:)); abs(y(:))]);
  close = isequal(size(x),size(y)) && all(abs(x(:)-y(:)) < tol*M);
end
