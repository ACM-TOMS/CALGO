% FIND_SM  Make the matrix Sm and its derivative
%
%    Sm = find_Sm(Scol, miss) returns Sm. SCol = [S0; S1; ...; S(n-1)].
%
%    [Sm, Smd] = find_Sm(Scol, miss, Scold) also finds the derivatives. Scold
%    are derivatives of Scol.

function [Sm, Smd] = find_Sm(Scol, miss, Scold)
  [r,n] = size(miss);
  DIFF = nargout > 1;
  %
  % CONSTRUCT OBSERVED INDICES
  miss = reshape(miss, r, n); %% for safeties sake
  obs = ~miss;
  [ko, ro, km, rm] = find_missing_info(miss);
  M = km(n+1);
  %
  % BUILD A COLUMN WITH 2n-1 MATRICES, [S(n-1)';...;S1';S0;S1;...S(n-1)]
  SU = reshape(Scol,r,n,r);
  SU = flipdim(permute(SU, [3,2,1]), 2);
  SS = [reshape(SU,r*n,r); Scol(r+1:end,:)];  
  Sm = zeros(M, M);
  if DIFF
    nPar = size(Scold,3);
    SUd = reshape(Scold,r,n,r,nPar);
    SUd = flipdim(permute(SUd, [3,2,1,4]), 2);
    SSd = cat(1, reshape(SUd,r*n,r,nPar), Scold(r+1:end,:,:));
    Smd = zeros(M, M, nPar);
  end
  %
  % CONSTRUCT Sm
  Sm = zeros(M, M);
  ib = r*(n-1)+1; ie = r*(2*n-1);
  for t=1:n
    if rm(t)~=0
      K = km(t)+1 : km(t+1);
      I = ib:ie;
      J = miss(:,t);
      Sm(:,K) = SS(I(miss), J);
      if DIFF, Smd(:,K,:) = SSd(I(miss), J, :); end
    end
    ib = ib - r;
    ie = ie - r;
  end
%   Sm = SomSm(miss,:);
%   if DIFF, Smd = SomSmd(miss,:,:); end
end
