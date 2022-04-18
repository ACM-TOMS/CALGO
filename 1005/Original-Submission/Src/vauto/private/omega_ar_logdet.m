% OMEGA_AR_LOGDET  Calculate log of determinant of Omega in pure AR case
%
%  LD = OMEGA_AR_LOGDET(Luo, LSig, jpat) finds log(det(Omega)) in the pure
%  autoregressive case. Luo is Cholesky factor of Su(obs,obs) and LSig{jpat(i)}
%  is the Cholesky factor of the i-th diagonal block in the lower right part of
%  Omega.
%
%  [LD, LDD] = OMEGA_AR_LOGDET(Luo, LSig, jpat, Luod, LSigd) finds also the
%  gradient of LD.

function [ld, ldd] = omega_ar_logdet(Luo, LSig, jpat, Luod, LSigd)
  npat = length(LSig);
  count = histc(jpat,1:npat);  % count of each pattern
  if nargout == 1  % ONLY FUNCTION VALUES
    ld = 2*sum(log(diag(Luo)));
    for j = 1:npat
      ld = ld + 2*count(j)*sum(log(diag(LSig{j})));
    end
  else  % FIND ALSO DERIVATIVES
    [ld, ldd] = logdet_L(Luo, Luod);
    ld = 2*ld; ldd = 2*ldd;
    for j = 1:npat
      [ldt, lddt] = logdet_L(LSig{j}, LSigd{j});
      ld = ld + 2*count(j)*ldt;
      ldd = ldd + 2*count(j)*lddt;
    end
  end
end
