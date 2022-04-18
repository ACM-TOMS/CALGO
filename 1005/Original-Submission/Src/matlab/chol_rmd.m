% CHOL_RMD  Reverse mode derivative of Cholesky factorization
%
%   Aa = CHOL_RMD(A,L,Aa,La) adds to Aa (the adjoint matrix of A) the
%   contribution from the Cholesky factorization of A (where f is a final output
%   result). L should be the Cholesky factor of A (computed earlier with L =
%   chol(A,'lower')) and La should be its (lower triangular) adjoint matrix.
%
%   For variety, the method given here is different from the BLAS-RMD Fortran
%   version, in spotrf_rmd.f90. The algorithm below is based on "compact"
%   Cholesky factorization, or the Cholesky-Crout algorithm, illustrated below,
%   whereas the Fortran version is based on the "recursive algorithm" (see
%   comments in the .f90 file). Both algorithms have the same time complexity,
%   but an advantage of the recursive-algorithm-based version is that it only
%   needs L and Aa as input, whereas the version below needs four input
%   arguments.
%
%   CHOLESKY FACTORIZATION CALCULATIONS (a, b, d AND s ARE IN COLUMN k):
%   ----A-----     ----L-----     ---step-k---
%   .              .              c = a - x'*x
%   :  .           :  .           d = sqrt(c)
%   x'    a        :     d        r = b - Y*x
%   Y     b  .     :     s  .     s = r/d

function Aa = chol_rmd(A, L, Aa, La)
  n = size(A,1);
  for k=n:-1:1
    d = L(k,k);
    s = L(k+1:n,k);
    x = L(k,1:k-1)';
    Y = L(k+1:n,1:k-1);
    %a = A(k,k);
    %b = A(k+1:n,k);   
    da = La(k,k);
    sa = La(k+1:n,k);
    aa = Aa(k,k);
    ba = Aa(k+1:n,k);
    xa = La(k,1:k-1)';
    Ya = La(k+1:n,1:k-1);
    
    % The actual computation:
    ra = sa/d;    
    da = da - s'*ra;
    ba = ba + ra;
    Ya = Ya - ra*x';
    xa = xa - Y'*ra;
    ca = da/(2*d);
    aa = aa + ca;
    xa = xa - 2*ca*x;
    
    Aa(k,k) = aa;
    Aa(k+1:n,k) = ba;
    La(k,1:k-1) = xa';
    La(k+1:n,1:k-1) = Ya;
  end
end
