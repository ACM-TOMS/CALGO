% CHOLF  Cholesky factorization with error message
%
%   L = CHOLF(A) finds the Cholesky factorization A = L*L'. Only the lower
%   triangle of A is accessed. If a non-positve definite A is encountered, an
%   error message that may be approriate for the VAR/VARMA functions is issued.

function L = cholf(A)
  [L,p] = chol(A');
  if p > 0    
    error('Cholesky factorization failed. Model may be non-stationary');
  else
    L = L';
  end
end
 