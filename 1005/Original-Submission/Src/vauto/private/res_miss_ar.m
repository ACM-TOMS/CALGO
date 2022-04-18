% RES_MISS_AR  Estimate residuals and missing values for AR-model
%
%   [EPS] = RES_MISS_AR(A,Sig,Lu,LSig,z) returns the maximum likelihood estimate
%   of the residuals of the VAR model:
%
%            x(t) - mu = A1·(x(t-1) - mu) + ... + Ap·(x(t-p) - mu) + eps(t)
%
%   where eps(t) is N(mu, Sig), mu is an r-vector, Sig is r×r and the A's are
%   the coefficients of the model (p>0). The input parameters are the variables
%   with the same names in var_llc.
%
%   [EPS] = RES_MISS_AR(A,Sig,Lu,LSig,woh,LR,LQ,K,u,v,Laomh,Vhat,Sm,miss) is
%   used for the missing value case. [EPS, XM] = FIND_AR_RES_MISS(...) returns
%   also the maximum likelihood estimate of missing values.

function [eps,xm] = res_miss_ar(A,Sig,Lu,LSig,z,LR,LQ,K,u,v,Laomh,Vhat,Sm,miss);
  LTT.LT = true; LTT.TRANSA = true;
  C = find_CGW(A, [], Sig);
  if nargin == 5
    r = size(Sig, 1);
    p = length(A)/r;
    n = length(z)/r;
    jpat = ones(1,n-p);
    y = omega_ar_trisolve(Lu, LSig, z, jpat, 0:r:n*r, 'T');
    y = lambda_T_multiply(A, y, r, n);
    miss = zeros(r, n);
    if nargout>1, xm = []; end
  else
    [r, n] = size(miss);
    p = size(A,2)/r;
    ko = find_missing_info(miss);
    [pat, ipat, jpat] = unique(miss(:, p+1:end)', 'rows');
    ko = find_missing_info(miss);
    v1 = linsolve(LQ, v, LTT);
    v2 = linsolve(LR, u + K*v1, LTT);
    % Following is because Vhat and Laomh are not stored in full if there are
    % no missing values at the end of the time series.
    k1 = size(Vhat,1);
    k2 = size(Laomh,1);
    y = z;
    y(1:k1) = y(1:k1) + Vhat*(v1 - v2);
    y(1:k2) = y(1:k2) + Laomh*(Sm*v2);
    y = omega_ar_trisolve(Lu, LSig, y, jpat, ko, 'T');
    y = lambda_T_multiply(A, y, r, n, miss);
    if nargout>1, xm = Sm*v2; end
  end
  eps = C_multiply(C, y, n, miss);
  eps = reshape(eps, r, n);
end
