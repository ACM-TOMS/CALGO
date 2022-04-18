function [a, b] = MDB_domain(MP)

% Computation of end points of the domain
%
% INPUT
%   MP    : MDB-spline multi-patch
%
% OUTPUT
%   a     : left end point
%   b     : right end point

a = MP.P(1).U(1);
b = a + sum(arrayfun(@(P) P.U(end) - P.U(1), MP.P));

end
