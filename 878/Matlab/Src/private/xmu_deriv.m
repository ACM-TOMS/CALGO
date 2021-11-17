%XMU_DERIV  Derivatives of observed components of x minus mu
%
%  XMUD = XMU_DERIV(NPAR, MISS) returns derivatives of the observed components
%  of x - repmat(mu,n,1) w.r.t. all NPAR parameters, A1...Ap, B1...Bq, Sig, mu
%  (with matrices taken column by column using only the lower triangle of Sig).
%  NPAR should be (p+q+1)·r^2 + r·(r+3)/2 and MISS(i,t) should be true if the
%  value X(i,t) is missing. Only the derivatives w.r.t. mu are nonzero. 

function xmud = xmu_deriv(nPar, miss)
  [r, n] = size(miss);
  [ko, ro] = find_missing_info(miss);
  obs = ~miss;
  nobs = sum(obs(:));
  xmud = zeros(nobs, 1, nPar);
  for t=1:n
    J = find(~miss(:,t));
    for i=1:ro(t)
      xmud(ko(t)+i, 1, nPar-r+J(i)) = -1;
    end
  end
end
