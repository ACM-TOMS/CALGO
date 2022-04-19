% FIND_ACOL  Create a block column of A-matrices and optionally derivatives
%
%   Acol = FIND_ACOL(A, r) returns the matrix
%
%                         A1
%                         :
%                         Ap
%
%   where I is the r×r identity matrix. 
%
%   [Acol, Acold] = FIND_ACOL(A, r, nPar) returns also the derivatives of Acol
%   w.r.t. the A matrices (taken column by column) and some additional nPar -
%   r^2·p parameters (normally the columns of moving average parameter matrices
%   Bj and a covariance matrix Sig).

function [Acol, Acold] = find_Acol(A, r, nPar)
  DIFF = nargout > 1;
  Ac = makecell(A);
  p = length(Ac);
  if p==0
    Acol = zeros(0, r); 
  else
    Acol= cat(1,Ac{:});
  end
  if DIFF
    Acold = zeros(r*p, r, nPar);
    k = 1;
    for j = 1:p
      for c = 1:r
        for l = 1:r
          Acold(l+r*(j-1),c,k) = 1;
          k = k+1;
        end
      end
    end
  end
end
