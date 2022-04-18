function [y,uflow,oflow]=quad_inerfc(N,n,x)
%QUAD_INERFC Quadrature approximation of the repeated 
%   integral f_n(x) of the coerror function.
%QUAD_INERFC(N,n,X) computes the N-point Gaussian quadrature
%   approximation y of f_n(X) and outputs uflow=1 if underflow
%   occurs, and oflow=1 if overflow occurs.

global ab ab0
uflow=0; oflow=0;
uv=gauss(N,ab0); u=uv(:,1); v=uv(:,2);
if x>=0
  y0=sum(v.*u.^n.*exp(-2*x*u));
  y=2*exp(-x^2)*y0/(sqrt(pi)*factorial(n));
  if isinf(y0), oflow=1; end
  if y==0, uflow=1; end
else
  y0=sum(v.*u.^n.*exp(-(x^2+2*x*u)));
  y=2*y0/(sqrt(pi)*factorial(n));
  if isinf(y0), oflow=1; end
  if y==0, uflow=1; end
end
