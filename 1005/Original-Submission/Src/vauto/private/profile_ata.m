% PROFILE_ATA  Multiply A transposed with A where A is a profile-sparse matrix
%
%   C = PROFILE_ATA(A, ko, km, p, g) sets C to A'·A making use of known sparsity
%   in A which is such that A = Af(obs,miss) and the (i,j)-block of Af is known
%   to be zero for i > max(p,j+g). This function is used to form Vhat'·Vhat
%   (with g=q) and Laomh'·Laomh (with g=p) efficiently. The vectors ko and km
%   specify the sparsity structure, ko(t) is the number of observed and km(t)
%   the number of missing values before time t.
%
%   [C,Cd] = PROFILE_ATA(A, ko, km, p, g, Ad) finds also the derivative of C.

function [C,Cd] = profile_ata(A, ko, km, p, g, Ad);
  n = length(ko)-1;
  DIFF = nargout>1;
  mA = size(A,1);
  M = size(A,2);
  C = zeros(M,M);
  if DIFF, nPar = size(Ad,3); Cd = zeros(M,M,nPar); end
  blk = ceil(25*n/(M+1));
  for t1=1:blk:n
    t = min(n,t1+blk-1);
    if km(t+1)>km(t1)
      tmax = min(n,max(p,t+g));
      I = 1 : min(mA, ko(tmax+1));
      K = km(t1)+1:km(t+1);
      L = 1 : km(t+1);
      J = 1 : km(t1);
      C(K,L) = C(K,L) + A(I,K)'*A(I,L);
      if DIFF
        if ~isempty(J)
          Cd(K,J,:) = Cd(K,J,:) + atb_deriv(A(I,K),Ad(I,K,:), A(I,J),Ad(I,J,:));
        end
        Cd(K,K,:) = Cd(K,K,:) + ata_deriv(A(I,K), Ad(I,K,:));
      end
    end
  end
  % Make upper triangular part:
  C = tril(C)+tril(C,-1)';
  if DIFF
    for k=1:size(Cd,3), Cd(:,:,k) = tril(Cd(:,:,k))+tril(Cd(:,:,k),-1)'; end
  end
end
