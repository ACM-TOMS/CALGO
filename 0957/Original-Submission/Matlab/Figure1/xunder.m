function x=xunder(n)
%XUNDER Upper boundary for a natural domain.
%XUNDER(N) determines, for given N, the value X=X^*(N), and
%   an approximation for it, above which F_N(X) underflows.

global ab ab0
N=90; tol=.5e-5;
for x=10:40
  [y,uflow]=quad_inerfc(N,n,x);
  if uflow==1
    a=x-1; b=x;
    break
  end
end
%[a b]
%[y,uflow]=quad_inerfc(N,n,a); uflow
%[y,uflow]=quad_inerfc(N,n,b); uflow
ktol=ceil(log((b-a)/tol)/log(2));
for k=1:ktol
%  [a b]
%  [y,uflow]=quad_inerfc(N,n,a); uflow
%  [y,uflow]=quad_inerfc(N,n,b); uflow
  x=(a+b)/2;
  [y,uflow]=quad_inerfc(N,n,x);
  if uflow==0, a=x; else b=x; end
end
