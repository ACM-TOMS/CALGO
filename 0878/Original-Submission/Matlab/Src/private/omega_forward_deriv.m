%OMEGA_FORWARD_DERIV  Derivative of solution to L·X = Y or Lo·X = Y
%
%  Xd = OMEGA_FORWARD_DERIV(Lu, Ll, Lud, Lld, X, Yd, p, q, ko) calculates the
%  derivative of the solution X to L·X = Y where L is the Cholesky factor of
%  Omega, stored in Lu and Ll (as returned by omega_factor), Lud and Lld are the
%  derivatives of Lu and Ll (calculated with omega_factor_deriv) and Yd is the
%  derivative of Y. ko(t) is the number of observed values before time t (for
%  the complete data case it is 0:r:r*n).
%
%  For the missing value case, Lu, Ll, Lud and Lld store the matrix Lo and its
%  derivatives.
%
%  As in OMEGA_FORWARD the calculation may be speeded when Y (and X) have a zero
%  block with ko(tmin) rows at top. If Y = [O; Y1], X = [O; X1] and Y1d is the
%  derivative of Y1, use the call X1d = OMEGA_FORWARD_DERIV(Lu, Ll, Lud, Lld,
%  X1, Y1d, p, q, ko, tmin). Also, to use only the leading tmax × tmax blocks of
%  L and its derivatives, call: OMEGA_FORWARD_DERIV(..., tmin, tmax).
%
%  METHOD: Differentiating through L·X = Y with respect to a parameter one
%  obtains L·Xd = Yd - Ld·X, so Xd may be obtained with forward substitution (by
%  calling omega_forward).
%    
function Xd = omega_forward_deriv(Lu, Ll, Lud, Lld, X, Yd, p, q, ko,tmin,tmax,T)
  if nargin<10, tmin=1; end
  if nargin<11, tmax=length(ko)-1; end
  n = length(ko) - 1;
  h = max(p,q);
  ro = diff(ko);
  e = ko(h+1);  % order of Lu
  Xd = Yd;
  nPar = size(Xd,3);
  for i=1:nPar
    % Find Yd-Ld·X and store in Xd:
    j = ko(tmin) + 1;
    m = ko(h+1) - j + 1;
    Xd(1:m,:,i) = Xd(1:m,:,i) - tril(Lud(j:end,j:end,i))*X(1:m,:);
    for t=max(h+1,tmin):tmax
      t1 = max(t-q,tmin);
      J = ko(t1)+1 : ko(t);
      K = ko(t)+1 : ko(t+1);
      JX = J - j + 1;
      KX = K - j + 1;
      JL = J - ko(t-q);
      KL = K - ko(t-q);
      Xd(KX,:,i) = Xd(KX,:,i)-Lld(K-e,JL,i)*X(JX,:)-tril(Lld(K-e,KL,i))*X(KX,:);
    end
    Xd(:,:,i) = omega_forward(Lu, Ll, Xd(:,:,i), p, q, ko, tmin, tmax);
  end
end
