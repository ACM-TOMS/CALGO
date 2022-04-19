function Hbar = MDB_nullspace(ll)

% Computation of left null-space of a column of matrix L
% (auxiliary function for MDB_extraction)
%
% INPUT
%   ll    : a column of L
%   
% OUTPUT
%   Hbar  : null-space matrix of ll

q = length(ll); 
i1 = find(ll, 1, 'first'); 
i2 = find(ll, 1, 'last');
dd = zeros(q-1, 2);
dd(1:i1, 1) = 1;
for j = i1:i2-2
   dd(j, 2) = -ll(j) / ll(j+1) * dd(j, 1);
   dd(j+1, 1) = 1 - dd(j, 2);
end
dd(i2-1:q-1, 2) = 1;
Hbar = spdiags(dd, [0 1], q-1, q);

end
