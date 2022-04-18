%OMEGA_LOGDET  Calculate log of determinant of Omega
%
%  LD = OMEGA_LOGDET(Lu, Ll, p, q, ko) returns log(det(Omega)) from the Cholesky
%  factorization of Omega given by omega_factor.
%
%  [LD, LDD] = OMEGA_LOGDET(Lu, Ll, p, q, ko, Lud, Lld) returns also the
%  derivative of log(det(Omega)) in LDD. Lud and Lld are from
%  omega_factor_deriv.
%
%  In the missing value case, Lu, Ll, Lud and Lld store the matrix L_o and its
%  derivative.

function [ld,ldd] = omega_logdet(Lu, Ll, p, q, ko, Lud, Lld)
  n = length(ko) - 1;
  if nargin <=5  % no derivatives
    nPar = 0;
    ld = 2*sum(log(diag(Lu)));
  else
    nPar = size(Lld,3);
    [ld, ldd] = logdet_L(Lu, Lud);
    ld = 2*ld; ldd = 2*ldd;
  end
  e = size(Lu,1);
  h = max(p,q);
  ro =  diff(ko);
  for t = h+1:n
    J = ko(t)-e + (1:ro(t));
    K = ko(t)-ko(t-q) + (1:ro(t));
    if nPar==0
      ld = ld + 2*logdet_L(Ll(J,K));
    else
      [ldt, lddt] = logdet_L(Ll(J,K), Lld(J,K,:));
      ld = ld + 2*ldt;
      ldd = ldd + 2*lddt;
    end
  end
end
