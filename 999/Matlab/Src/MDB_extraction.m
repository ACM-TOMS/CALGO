function H = MDB_extraction(MP, rr)

% Computation of multi-degree spline extraction matrix
%
% INPUT
%   MP    : MDB-spline multi-patch
%   rr    : MDB-spline smoothness vector (optional)
%
% OUTPUT
%   H     : extraction matrix

if nargin < 2
   rr = 0;
end
m = length(MP.P) - 1;
if length(rr) == 1
   rr = repmat(rr, 1, m);
end
H = speye(MP.mu(end));
for i = 1:m
   r = min([rr(i), MP.P(i).p, MP.P(i+1).p]);
   K = sparse([B_diffend_all(MP.P(i), r, false);
              -B_diffend_all(MP.P(i+1), r, true)]);
   L = H(:, MP.mu(i)+1:MP.mu(i+2)) * K;
   for j = 0:r
      Hbar = MDB_nullspace(L(:, j+1));
      H = Hbar * H;
      L = Hbar * L;
   end
end

end
