% RES_MISS  Estimate residuals and missing values
%
%   EPS = RES_MISS(A, C, Lu, Ll, w) returns the maximum likelihood
%   estimate of the residuals of the VARMA model:
%
%            x(t) - mu = A1·(x(t-1) - mu) + ... + Ap·(x(t-p) - mu) + y(t)
%   with
%            y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q),
%
%   where eps(t) is N(mu, Sig), mu is an r-vector, Sig is r×r and the A's and
%   B's are the coefficients of the model. Lu, Ll, C and w are as in varma_llc.
%
%   [EPS, XM] = RES_MISS(A,C,Luo,Llo,wohat,LR,LQ,K,u,v,Laomh,Vhat,Sm,miss)
%   is used for the missing value case. The paramteres are as in varma_llm, and
%   in addition to the EPS estimate, the maximum likelihood estimate of the
%   missing values is returned in XM.

function [eps,xm] = res_miss(varargin);
  if nargin==5
    [A, C, Lu, Ll, w] = deal(varargin{:});
    r = size(C{1}, 1);
    p = length(A)/r;
    q = length(C) - 1;
    n = length(w)/r;
    what = omega_forward(Lu, Ll, w, p, q, 0:r:r*n);
    y = omega_back_sub(Lu, Ll, what, p, q, 0:r:r*n);
    y = lambda_T_multiply(A, y, r, n);
    eps = C_multiply(C, y, n);
    xm = [];
  else
    [A, C, Luo, Llo, wohat, LR,LQ,K,u,v,Laomh,Vhat,Sm,miss] = deal(varargin{:});
    [r, n] = size(miss);
    p = length(A)/r;
    q = length(C) - 1;
    ko = find_missing_info(miss);
    LTT.LT = true; LTT.TRANSA = true;
    v1 = linsolve(LQ, v, LTT);
    v2 = linsolve(LR, u + K*v1, LTT);
    % Following is because Vhat and Laomh are not stored in full if there are
    % no missing values at the end of the time series.
    k1 = size(Vhat,1);
    k2 = size(Laomh,1);
    y = wohat;
    y(1:k1) = y(1:k1) + Vhat*(v1 - v2);
    y(1:k2) = y(1:k2) + Laomh*(Sm*v2);
    y = omega_forward(Luo, Llo, y, p, q, ko);
    y = lambda_T_multiply(A, y, r, n, miss);
    eps = C_multiply(C, y, n, miss);
    xm = Sm*v2;
  end
  eps = reshape(eps, r, n);
end
