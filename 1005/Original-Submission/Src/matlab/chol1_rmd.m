% Matlab version of the Cholesky adjoint that matches the Fortran version
function La = chol1_rmd(L, La)
  n = size(L,1);
  K = 2:n;
  if n > 1
    La(K,K) = chol1_rmd(L(K,K), La(K,K));
    La(K,1) = (La(K,1) - (La(K,K) + La(K,K)')*L(K,1))/L(1,1);
  end
  La(1,1) = (La(1,1) - L(K,1)'*La(K,1))/(2*L(1,1));
end
