%OMEGA_FACTOR_DERIV  Derivative of Cholesky factorization of Omega
%
%  [Lud,Lld] = OMEGA_FACTOR_DERIV(Sud, Olowd, Lu, Ll, p, q, ko) calculates the
%  derivative of the Cholesky factorization of Omega. On entry Sud should be an
%  r·h×r·h×nPar array with the derivatives of Su and Olowd should be an
%  r·(n-h)×r·(q+1) array with the derivatives of Olow where Su and Olow are the
%  upper left and lower partitions of Omega as returned by omega_build.
%  Moreover, Lu and Ll should be as returned by omega_factor and p and q are the
%  dimensions of the model. ko is an (n+1)-vector with ko(t) = number of
%  observed values before time t. On exit Lud and Lld return the derivatives of
%  Lu and Ll.
%    
%  In the missing value case, Sud and Olowd are from omega_remove_miss.

function [Lud, Lld] = omega_factor_deriv(Sud, Olowd, Lu, Ll, p, q, ko)
  n = length(ko) - 1;
  h = max(p, q);
  Lud = chol_deriv(Lu, Sud);  % Derivative of upper left partition:
  %
  % Find derivative of lower part:
  e = size(Lu, 1);
  Lld = zeros(size(Olowd));
  for t = h+1 : n  % loop over block-lines in Olow
    K = ko(t)+1 : ko(t+1);    
    KL = K - ko(t-q);
    JL = 1 : ko(t)-ko(t-q);
    U = Ll(K-e, JL);
    Yd = permute(Olowd(K-e, JL, :), [2,1,3]);
    Ud = omega_forward_deriv(Lu, Ll, Lud, Lld, U', Yd, p, q, ko, t-q, t-1);
    Ud = permute(Ud, [2,1,3]);
    Lld(K-e, JL, :) = Ud;
    UUdT = aat_deriv(U, Ud);
    Lld(K-e, KL, :) = chol_deriv(Ll(K-e, KL), Olowd(K-e, KL, :) - UUdT);
  end
end    
