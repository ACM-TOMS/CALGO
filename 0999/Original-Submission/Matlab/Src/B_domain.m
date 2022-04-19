function [a, b] = B_domain(P)

% Computation of end points of the domain
%
% INPUT
%   P     : B-spline patch
%
% OUTPUT
%   a     : left end point
%   b     : right end point

a = P.U(1);
b = P.U(end);

end
