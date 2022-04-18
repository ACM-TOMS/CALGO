% FIND_LAMBDA_OM  Find the matrix Lambda_om
%
%   Lamom = FIND_LAMBDA_OM(Acol, miss) returns the matrix consisting of observed
%   rows and missing columns of Lambda i.e. Lam(obs,miss) where Lam is the
%   matrix shown in the help of lambda_multiply. Acol should be [A1; ... Ap],
%   provided by find_Acol.
%   
%   [Lamom,Lamomd] = FIND_LAMBDA_OM(Acol, miss, Acold) finds also the derivative
%   (Acold contains the derivatives of Acol).

function [Laom,Laomd] = find_lambda_om(Acol, miss, Acold)
  DIFF = nargout > 1;
  [r,n] = size(miss);
  obs = ~miss;
  [ko,ro,km,rm] = find_missing_info(miss);
  N = ko(n+1);
  M = km(n+1);
  p = size(Acol,1)/r;
  mLaom = findmaxnz(rm, ko, p, n);
  Laom = zeros(mLaom, M);
  if DIFF
    nPar = size(Acold,3);
    Laomd = zeros(mLaom, M, nPar);
  end
  for t=1:n
    if rm(t)~=0
      tmax = min(n, t+p);
      mY = ko(tmax+1);
      I = obs(:,t+1:tmax);
      J = miss(:,t);
      K = km(t)+1:km(t+1);
      Laom(ko(t+1 )+1:mY, K) = -Acol(I,J);
      if t<p, Laom(ko(t+1)+1:ko(p+1), K) = 0; end
      if DIFF
        Laomd(ko(t+1)+1:mY, K, :) = -Acold(I,J,:);
        if t<p, Laomd(ko(t+1)+1:ko(p+1), K, :) = 0; end
      end
    end
  end
end

function maxnz = findmaxnz(rm, ko, p, n)
  % return maximum non-zero row-number of Lambda_om
  tmax = find(rm,1,'last');  % last missing time
  trow = min(n, tmax + p);   % time of last nonzero row
  maxnz = ko(trow + 1);      % number of observations until that time
end
