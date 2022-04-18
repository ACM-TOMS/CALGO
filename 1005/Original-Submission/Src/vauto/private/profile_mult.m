% PROFILE_MULT  Multiply with a profile-sparse matrix
%
%   C = PROFILE_MULT(A, B, ko, km, p, g, h) sets C to A'*B making use of known
%   sparsity in A and B which is such that A = Af(obs,miss) and the (i,j)-block
%   of Af is known to be zero for i > max(p,j+g) and B = Bf(obs,miss) with the
%   (i,j)-block of Bf zero for i > max(p,j+h). When B is full set h = n. This
%   function is used to multiply efficiently with Vhat' (g=q) and Laomh' (g=p).
%   The matrix P = Laomh'·Vhat is formed with g=p, h=q. ko(t) is the number of
%   observed and km(t) the number of missing values before time t.
%
%   [C,Cd] = PROFILE_MULT(A,B, ko, km, p, g, h, Ad,Bd) finds also the derivative
%   of C.

function [C,Cd] = profile_mult(A, B, ko, km, p, g, h, Ad, Bd);
  n = length(ko)-1;
  DIFF = nargout>1;
  [k1,l] = size(A); [k2,M] = size(B);
  m = min(k1, k2);
  C = zeros(l,M);
  if DIFF, nPar = size(Ad,3); Cd = zeros(l,M,nPar); end
  blk = ceil(25*n/(M+1));
  if h<n
    tA1 = 1;
    tB1 = 1;
    while tA1<=n && tB1<=n
      tA = min(n,tA1 + blk - 1);
      tB = min(n,tB1 + blk - 1);
      kA = min(m, ko(min(n, max(p,tA+g)) + 1));
      kB = min(m, ko(min(n, max(p,tB+h)) + 1));
      if kA < kB
        if km(tA+1) > km(tA1)
          I = km(tA1)+1:km(tA+1);
          J = 1:kA;
          L = km(tB1)+1:km(n+1);
          C(I,L) = C(I,L) + A(J,I)'*B(J,L);
          if DIFF
            Cd(I,L,:) = Cd(I,L,:)+atb_deriv(A(J,I),Ad(J,I,:),B(J,L),Bd(J,L,:));
          end
        end
        tA1 = tA1+blk;
      else
        if km(tB+1) > km(tB1)
          I = km(tB1)+1:km(tB+1);
          J = 1:kB;
          L = km(tA1)+1:km(n+1);
          C(L,I) = C(L,I) + A(J,L)'*B(J,I);
          if DIFF
            Cd(L,I,:) = Cd(L,I,:)+atb_deriv(A(J,L),Ad(J,L,:),B(J,I),Bd(J,I,:));
          end
        end
        tB1 = tB1+blk;
      end
    end
  else
    for t1=1:blk:n
      t = min(n,t1+blk-1);
      if km(t+1) > km(t1)
        tmax = min(n,max(p,t+g));
        I = 1 : min(m, ko(tmax+1));
        K = km(t1)+1:km(t+1);
        C(K,:) = C(K,:) + A(I,K)'*B(I,:);
        if DIFF
          XX = atb_deriv(A(I,K),Ad(I,K,:), B(I,:),Bd(I,:,:));
          Cd(K,:,:) = Cd(K,:,:) + XX;
        end
      end
    end
  end
end
