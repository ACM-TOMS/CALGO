% FIND_V  Find V-matrix (and its derivative) for missing value likelihood
%
%   V = FIND_V(G, Gcol, Scol, miss) finds the matrix V as the observed rows
%   and missing columns of the matrix Lam·S which is given by:
%
%        S0 S1'...        Sn-1'      S0 S1'...     Sn-1'
%        S1 S0               :       S1 S0            : 
%         :        ...      S1'       :       ...    S1'
%        Sp-1...            S0       Sp-1...         S0 
%        ----------------------  or  -------------------
%              Gq ...    Gp-n+1      Gp ...       Gp-n+1
%                 Gq         :        :               :
%                    ..              Gq       
%                      Gq...G0          ..
%                                         Gq...      G0
%
%   (depending on whether p or q is bigger). G, Gneg and Scol should be
%   {G0...Gq} (from find_CGW), [G(p-n+1);...;G(-1)] (from find_Gneg) and
%   [S0;...; S(n-1)] (from vyw_solve and S_extend). miss(i,t) is false if X(i,t)
%   is missing.
%
%   [V, Vd] = FIND_V(G, Gneg, Scol, miss, Gd, Gnegd, Scold) finds also the
%   derivative of V. Gd and Scold should be (cell) arrays with the derivatives
%   of G and Scol. 

function [V,Vd] = find_V(G, Gneg, Scol, miss, Gd, Gnegd, Scold)
  DIFF = nargout > 1;
  [r,n] = size(miss);
  obs = ~miss;
  [ko,ro,km,rm] = find_missing_info(miss);
  N = ko(n+1);
  M = km(n+1);
  p = n - size(Gneg,1)/r - 1;
  q = length(G) - 1;
  mV = findmaxnz(rm, ko, p, q, n);
  V = zeros(mV, M);
  if DIFF
    nPar = size(Scold,3);
    Vd = zeros(mV, M, nPar);
  end
  % Make SS = [S(n-1)';...; S2'; S1'; S0;... S(p-1)]:
  SU = reshape(Scol,r,n,r); 
  SU = flipdim(permute(SU, [3,2,1]), 2);
  SS = [reshape(SU,r*n,r); Scol(r+1:r*p,:)]; 
  GG = cat(1,Gneg,G{:});
  if DIFF
    SUd = reshape(Scold,r,n,r,nPar);
    SUd = flipdim(permute(SUd, [3,2,1,4]), 2);
    SSd = cat(1, reshape(SUd,r*n,r,nPar), Scold(r+1:r*p,:,:));
    GGd = cat(1,Gnegd,Gd{:});
  end
  for t=1:n
    if rm(t)~=0
      LamS = [
        SS((n-t)*r+1 : (n-t+p)*r, miss(:,t));
        GG((n-t)*r+1 : end, miss(:,t))];
      tmax = min(n,max(p,t+q));
      mY = ko(tmax+1);
      K = km(t)+1:km(t+1);
      V(1:mY, K) = LamS(obs(:,1:tmax), :);
      if DIFF
        LamSd = cat(1, ...
          SSd((n-t)*r+1 : (n-t+p)*r, miss(:,t), :), ...
          GGd((n-t)*r+1 : end, miss(:,t), :));
        Vd(1:mY, K, :) = LamSd(obs(:,1:tmax), :, :);
      end
    end
  end
end

function maxnz = findmaxnz(rm, ko, p, q, n)
  % return maximum non-zero row-number of V
  tmax = find(rm,1,'last');       % last missing time
  if isempty(tmax), tmax=0; end;  % in case nothing is missing
  trow = min(n, max(p, tmax+q));  % time of last nonzero row
  maxnz = ko(trow + 1);           % number of observations until that time
end
