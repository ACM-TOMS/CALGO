function gg = MDB_greville(MP, H)

% Computation of multi-degree Greville points
%
% INPUT
%   MP    : MDB-spline multi-patch
%   H     : extraction matrix
%
% OUTPUT
%   gg    : vector of Greville points

x = 0;
m = length(MP.P) - 1;
dd = zeros(1, MP.mu(end));
for i = 1:m
   dd(MP.mu(i)+1:MP.mu(i+1)) = x + B_greville(MP.P(i));
   x = x + MP.P(i).U(end) - MP.P(i+1).U(1);
end
dd(MP.mu(m+1)+1:MP.mu(end)) = x + B_greville(MP.P(end));
gg = dd / H;

end
