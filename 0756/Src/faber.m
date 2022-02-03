function F = faber(m,w,beta,z,c)
%FABER  Faber polynomial coefficients for polygonal regions.
%       FABER(M,W,BETA,Z,C) returns the coefficients of Faber
%       polynomials of degree <= M for the polygonal region described by
%       W and BETA.  A call to DEPARAM must be made first to obtain the
%       values of Z and C for the Schwarz-Christoffel exterior map.
%       FABER will return an upper triangular square matrix P of size
%       M+1 such that P(1:k,k) is the vector of coefficients for the
%       Faber polynomial of degree k-1.  Note that the leading (highest
%       degree) coefficient is always first.
%	
%	See also DEPARAM, FABERDEMO.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

%       This function follows somewhat closely the procedure outlined in
%       section 4 of Starke and Varga (Num. Math., 1993), except that no
%       symmetry of the polygon is assumed.

if nargin < 6
  qdat = scqdata(beta,8);
end
n = length(w);
gam = ones(n,m+1);			% coeffs of binomial expansion
Z = ones(n,m-1);			% powers of the z(j)
for k = 1:m
  gam(:,k+1) = -gam(:,k).*(beta-k+1)./(k*z);
  if k < m
    Z(:,k) = z.^k;
  end
end

% Compute the coeffs of the Laurent expansion of Psi
e1 = zeros(m+1,1);
e1(1) = 1;
C = -c*e1;
for j = 1:n
  C = toeplitz(gam(j,:).', e1')*C;
end
C = C(3:m+1)./(-(1:m-1)');
%c0 = (sum(w) + c*sum(1./z) - sum(Z)*C)/n;
x0 = 10^(-10/m);
c0 = demap(x0,w,beta,z,c,qdat) + c/x0 - x0.^(1:m-1)*C;
C = [c0;C];

% Use the Faber recurrence to compute polynomial coeffs
P = zeros(m+1,m+1);
P(1,1) = 1;				% poly coeffs, low order first
F = P;					% high order first (MATLAB style)
for k = 1:m
  P(1:k+1,k+1) = ([0;P(1:k,k)] - P(1:k+1,1:k)*[k*C(k);C(k-1:-1:1)])/(-c);
  F(1:k+1,k+1) = flipud(P(1:k+1,k+1));
end

% Normalize so that leading coefficients are real
F = F*diag(exp(-i*angle(F(1,:))));

