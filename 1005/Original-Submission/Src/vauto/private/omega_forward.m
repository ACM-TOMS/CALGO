% OMEGA_FORWARD  Forward substitution with Cholesky factor of Omega or Omega_o
%
%  X = OMEGA_FORWARD(Lu, Ll, Y, p, q, ko) solves L·X = Y where p, q and r are
%  the dimensions of the ARMA model and L is the Cholesky factor of Omega,
%  stored in Lu and Ll (as returned by omega_factor). ko(t) is the number of
%  observed values before time t. In the missing value case, Lu and Ll store the
%  matrix L_o (see omega_factor).
%
%  When Y = [O; Y1] where O is a zero matrix with ko(tmin) rows, the calculation
%  may be speeded using the call X1 = OMEGA_FORWARD(Lu, Ll, Y1, p, q, ko, tmin).
%  The X1 returned will contain rows ko(tmin)+1...ko(n+1) of X (the first
%  ko(tmin) rows being known to be zero). Moreover, to use only the first tmax
%  block rows of L, use X1 = OMEGA_FORWARD(Lu, Ll, Y1, p, q, ko, tmin, tmax),
%  with Y1 containing block-rows tmin to tmax. This call syntax is used in
%  omega_factor.

function X = omega_forward(Lu, Ll, Y, p, q, ko, tmin, tmax)
  LT.LT = true;
  if nargin<7, tmin=1; end
  if nargin<8, tmax=length(ko)-1; end
  h = max(p,q);
  ro = diff(ko);
  e = ko(h+1);  % order of Lu 
  j = ko(tmin) + 1;
  m = ko(h+1) - j + 1;
  X = Y;
  if m>0, X(1:m,:) = linsolve(Lu(j:end,j:end), Y(1:m,:), LT); end
  for t=max(h+1,tmin):tmax
    t1 = max(t-q,tmin);
    J = ko(t1)+1 : ko(t);
    K = ko(t)+1 : ko(t+1);
    JX = J - j + 1;
    KX = K - j + 1;
    JL = J - ko(t-q);
    KL = K - ko(t-q);
    X(KX,:) = X(KX,:) - Ll(K-e,JL)*X(JX,:);
    X(KX,:) = linsolve(Ll(K-e,KL), X(KX,:), LT);
  end
end
