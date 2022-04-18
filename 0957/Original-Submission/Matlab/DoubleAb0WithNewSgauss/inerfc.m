function [y,N]=inerfc(n,x)
%INERFC Repeated integral of the coerror function.
%INERFC(n,X) computes the function y=f_n(X) and outputs the
%   number N of Gaussian quadrature points used.

global ab0
if n<0 | n>150, error('n out of range'); end
if x<-17.9+.024*n | x>27-.084*n, error('x out of range'); end
if x>=0
  if n<=25, N=ceil(21.5+2.388*x);
  elseif n<=50, N=ceil(24+2.08*x);
  elseif n<=75, N=ceil(33+1.549*x);
  elseif n<=100, N=ceil(40+1.068*x);
  elseif n<=125, N=ceil(47+.4348*x);
  else N=ceil(54-.1852*x);
  end
else
  N=ceil(10+10.5*abs(x)+(.3223+.00747*abs(x))*n);
  if N>200, N=200; end
end
uv=gauss(N,ab0); u=uv(:,1); v=uv(:,2);
if x>=0
  y0=sum(v.*u.^n.*exp(-2*x*u));
  y=2*exp(-x^2)*y0/(sqrt(pi)*factorial(n));
else
  y0=sum(v.*u.^n.*exp(-(x^2+2*x*u)));
  y=2*y0/(sqrt(pi)*factorial(n));
end
