function e = end(x, k, n)
%END Last index in one mode of htensor.
%
%   Example
%   The expression X(end,:,:) calls END(X,1,3) to determine the value of
%   the last index in mode 1.
%
%   See also HTENSOR, HTENSOR/SUBSREF, END.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Check that n == ndims(x)
if( n ~= ndims(x) )
  error('Wrong number of indices: ndims(X) = %d.', ndims(x));
end

e = size(x, k);
