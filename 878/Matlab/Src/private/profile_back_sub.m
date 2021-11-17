% PROFILE_BACK_SUB  Omega back substitution for a profile-sparse matrix
%
%   X = PROFILE_BACK_SUB(Lu, Ll, Y, ko, km, p, q, g) solves Lo'·X = Y exploiting
%   known sparsity in Y, which is such that Y = Yf(obs,miss) and the (i,j)-block
%   of Yf is known to be zero for i > max(p,j+g). This function is used to form
%   Vhat (with g=q) and Laomh (with g=p). Lo is the ltl-factor of Omega_o (from
%   omega_ltl), stored in Lu and Ll, and ko(t) is the number of observed and
%   km(t) the number of missing values before time t..
%
%   [X,Xd] = PROFILE_BACK_SUB(Lu, Ll, Y, ko, km, p, q, g, Lud, Lld, Yd) finds
%   also the derivative of X.

%   NOTE: The routine is implemented using a blocked algorithm for speed

function [X,Xd] = profile_back_sub(Lu, Ll, Y, ko, km, p, q, g, Lud, Lld, Yd)
  DIFF = nargout > 1;
  n = length(ko) - 1;
  M = km(n+1);
  blk = ceil(25*n/(M+1));
  X = Y;
  if DIFF, Xd = Yd; end
  mY = size(Y,1);
  for t1=1:blk:n
    t = min(n,t1+blk-1);
    if km(t+1) > km(t1)
      tmax = min(n,max(p,t+g));
      I = 1 : min(mY, ko(tmax+1));
      K = km(t1)+1:km(t+1);
      if ~DIFF
        X(I,K) = omega_back_sub(Lu, Ll, Y(I,K), p, q, ko);
      else
        [X(I,K),Xd(I,K,:)] = ...
          omega_back_sub(Lu, Ll, Y(I,K), p, q, ko, Lud, Lld,Yd(I,K,:));
      end
    end
  end
end
