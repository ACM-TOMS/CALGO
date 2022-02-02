function ss = B_evaluation_spline(P, cc, xx)

% Evaluation of a spline in given points
%
% INPUT
%   P     : B-spline patch
%   cc    : vector of coefficients
%   xx    : vector of evaluation points
%
% OUTPUT
%   ss    : vector of spline evaluation values

M = B_evaluation_all(P, xx);
ss = reshape(cc, 1, []) * M;

end
