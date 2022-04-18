function ss = B_differentiation_spline(P, r, cc, xx)

% Differentiation of a spline in given points
%
% INPUT
%   P     : B-spline patch
%   r     : order of derivative
%   cc    : vector of coefficients
%   xx    : vector of evaluation points
%
% OUTPUT
%   ss    : vector of r-th order derivative spline values

M = B_differentiation_all(P, r, xx);
ss = reshape(cc, 1, []) * M;

end
