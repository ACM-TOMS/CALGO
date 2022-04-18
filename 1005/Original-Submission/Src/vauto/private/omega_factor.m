%OMEGA_FACTOR  Cholesky factorization of Omega
%
%  [Lu,Ll,info] = omega_factor(Su,Olow,p,q,ko) calculates the Cholesky
%  factorization Omega = L·L' of Omega which is stored in two parts, a full
%  upper left partition, Su, and a block-band lower partition, Olow, as returned
%  by omega_build. Omega is symmetric, only the lower triangle of Su is
%  populated, and Olow only stores diagonal and subdiagonal blocks. On exit, L =
%  [L1; L2] with L1 = [Lu 0] and L2 is stored in block-band-storage in Ll. Info
%  is 0 on success, otherwise the loop index resulting in a negative number
%  square root. P and q are the dimensions of the problem and ko is a vector
%  with ko(t) = number of observed values before time t.
%
%  In the complete data case ko should be 0:r:n*r. For missing values, Su and
%  Olow are the upper left and lower partitions of Omega_o = Omega with missing
%  rows and columns removed. In this case Lu and Ll return L_o, the Cholesky
%  factor of Omega_o.

function [Lu,Ll,info] = omega_factor(Su,Olow,p,q,ko)
  n = length(ko)-1;
  h = max(p,q);
  ro = diff(ko);
  Ll = zeros(size(Olow));
  [Lu,info] = chol(Su');  % upper left partition
  if info>0, return; end
  Lu = Lu';
  e = ko(h+1);     % order of Su
  for t = h+1 : n  % loop over block-lines in Olow
    K = ko(t)+1 : ko(t+1);    
    KL = K - ko(t-q);
    JL = 1 : ko(t)-ko(t-q);
    tmin = t-q; tmax = t-1;
    Ll(K-e, JL) = omega_forward(Lu, Ll, Olow(K-e,JL)', p, q, ko, tmin, tmax)';
    [Ltt, info] = chol(Olow(K-e,KL) - Ll(K-e,JL)*Ll(K-e,JL)');
    if info>0, info = info + ko(t); return; end
    Ll(K-e, KL) = Ltt';
  end
end
