% GAUSSRAT_RATIONAL Rational Gauss quadrature rule.
%
%    The call xw=GAUSS_RATIONAL(N,abmod) generates from the first N 
%    modified recurrence coefficients the N-point rational Gauss-type 
%    quadrature xw according to Theorem 3.25.
%
function xw=gauss_rational(N,abmod)
global Z M
xw=gauss(N,abmod);
if M==0, return, end
for n=1:N
  p(n)=prod((1+xw(n,1)*Z(1:M,1)).^Z(1:M,2));
end
xw(:,2)=xw(:,2).*p';

