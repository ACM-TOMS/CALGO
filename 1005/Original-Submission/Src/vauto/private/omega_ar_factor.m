% OMEGA_AR_FACTOR  Cholesky factorization of Omega in pure autoregressive case
%
%  [Lu, LSig, jpat] = OMEGA_AR_FACTOR(S, Sig, miss) determines the Cholesky
%  factorization of Omega(obs,obs) in the pure autregressive case. Omega is then
%
%    S0 S1'...      Sp-1'|
%    S1 S0             : |
%    S2 S1 S0            |             r·p
%     :        ...    S1'|
%    Sp-1             S0 |
%    --------------------+------------
%                        |Sig
%                        |   Sig       r·(n-p)
%                        |      ...
%                        |         Sig
%
%  Lu returns the Cholesky factor of the upper left partition and LSig{j}
%  returns the Cholesky factor of Sig(obs,obs) where obs are the indices of the
%  observed values for the j-th missing pattern occurring after time p.
%  jpat(i) returns the number of the missing pattern at time p+j. Thus
%  LSig{jpat(i)} is the Cholesky factor of the i-th block in the lower right
%  partition of Omega(obs,obs). S should be {S0...S(p-1)} and miss should be
%  r×n, true in missing locations.
%
%  [Lu, LSig, jpat, Lud, LSigd] = OMEGA_AR_FACTOR(S, Sig, miss, Sd, Sigd)
%  finds also the derivatives of Lu and LSig.
%
%  OMEGA_AR_FACTOR(S, Sig, n) and OMEGA_AR_FACTOR(S, Sig, n, Sd, Sigd) may be
%  used instead when there are no missing values.

function [Lu,LSig,jpat,Lud,LSigd] = omega_ar_factor(S,Sig,miss,Sd,Sigd)
  DIFF = nargout > 3; if DIFF, ascertain(nargin==5); end
  MISS = length(miss) > 1 && any(miss(:));
  p = max(0, length(S) - 1);
  Su = omega_build(S, {Sig}, {Sig}, p, p);
  r = size(Sig,1);
  if MISS
    n = size(miss, 2);
    miss1 = miss(:,1:p);
    miss2 = miss(:,p+1:end);
  else
    if length(miss)==1, n = miss; else n = size(miss,2); end
  end
  if DIFF
    Sud = omega_build_deriv(Sd, {Sigd}, {Sigd}, p, p);
    if MISS,
      [Su, Olow, Sud] = omega_remove_miss(Su, [], miss1, Sud, []);
    end
  elseif MISS
    Su = omega_remove_miss(Su, [], miss1);
  end
  Lu = cholf(Su);
  if DIFF, Lud = chol_deriv(Lu, Sud); end
  if MISS
    [pat, ipat, jpat] = unique(miss2', 'rows'); % the unique missing patterns
  else
    pat = false(1, r);
    jpat = ones(1, n - p);
  end
  npat = size(pat,1);
  for i=1:npat
    obs = ~pat(i,:);
    LSig{i} = cholf(Sig(obs,obs));
    if DIFF
      LSigd{i} = chol_deriv(LSig{i}, Sigd(obs, obs, :)); 
    end
  end
end
