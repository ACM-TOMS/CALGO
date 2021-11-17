function gg = B_greville(P)

% Computation of Greville points
%
% INPUT
%   P     : B-spline patch
%
% OUTPUT
%   gg    : vector of Greville points

gg = mean(P.U(hankel(2:P.p+1, P.p+1:P.p+P.n)), 1);

end
