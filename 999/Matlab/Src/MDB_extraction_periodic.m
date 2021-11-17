function H = MDB_extraction_periodic(MP, rr, rp)

% Computation of multi-degree spline extraction matrix with periodicity
%
% INPUT
%   MP    : MDB-spline multi-patch
%   rr    : MDB-spline smoothness vector (optional)
%   rp    : periodicity smoothness (optional)
%           the value should be less than half the dimension (floored) 
%           of the related non-periodic MDB-spline space
%
% OUTPUT
%   H     : extraction matrix

if nargin < 2
   rr = 0;
end
H = MDB_extraction(MP, rr);
if nargin > 2
   m = length(MP.P);
   r = min([rp, MP.P(m).p, MP.P(1).p]);
   if size(H, 1) >= 2*(r+1)
      Hper = circshift(H, r+1);   
      K = sparse([B_diffend_all(MP.P(m), r, false);
                 -B_diffend_all(MP.P(1), r, true)]);
      Lper = Hper(:, [MP.mu(m)+1:MP.mu(m+1) 1:MP.mu(2)]) * K;
      for j = 0:r
         Hbar = MDB_nullspace(Lper(:, j+1));
         Hper = Hbar * Hper;
         Lper = Hbar * Lper;
      end
      H = Hper;
   end
end

end
