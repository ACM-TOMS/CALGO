function P = B_patch(p, xx, kk)

% Construction of a B-spline patch of degree p with open knot vector
%
% INPUT
%   p     : B-spline degree
%   xx    : vector of break points
%   kk    : smoothness vector (optional)
%
% OUTPUT
%   P     : B-spline patch
%
% If kk is a scalar, smoothness kk is imposed at break point xx(i+1),
% if kk is a vector, smoothness kk(i) is imposed at break point xx(i+1),
% for i = 1:length(xx)-2

if nargin < 3
   kk = 0;
end
if length(kk) == 1
   kk = repmat(kk, 1, length(xx)-2);
end
U = repelem(xx, [p+1, p-kk, p+1]);
n = length(U) - p - 1;
P = struct('p', p, 'n', n, 'U', U);

end
