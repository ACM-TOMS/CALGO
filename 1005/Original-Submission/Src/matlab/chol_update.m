function [A,Aa] = chol_update(A,La)
  n = size(A,1);
  % "Recursive" Cholesky factorization
  for k = 1:n
    J = k+1:n;
    a = A(k,k);
    v = A(J,k);
    B = A(J,J);
    a = sqrt(a);
    v = v/a;
    C = B - v'*v;
    A(k,k) = a;
    A(k,J) = 0;
    A(J,k) = v;
    A(J,J) = tril(C);
  end
  Aa = chol_update_rmd(A,La);
end
