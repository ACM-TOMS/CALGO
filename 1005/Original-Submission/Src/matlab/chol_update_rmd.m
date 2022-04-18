function Aa = chol_update_rmd(A,Aa)
  % Version of chol_rmd that updates La with Aa, assumes Aa=0 on entry,
  % needs only L and La (which are in A and Aa on entry). This version
  % corresponds to dpotrf_rmd.f90.
  n = size(A,1);
  Aa = tril(Aa);
  for k = n:-1:1
    J  = k+1:n;
    Ca = Aa(J,J);
    a  = A(k,k);
    v  = A(J,k);      
    aa = Aa(k,k);
    va = Aa(J,k);
    %
    va = (va - Ca*v - Ca'*v)/a;
    aa = (aa - va'*v)/(2*a);
    %
    Aa(J,k) = va;
    Aa(k,k) = aa;
  end
end
