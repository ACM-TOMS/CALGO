% PROFILE_AR_TRISOLVE  Pure AR Omega triangular solve for profile-sparse matrix
%
%   X = PROFILE_AR_TRISOLVE(Lu, LSig, Y, JPAT, ko, km, p, g, 'N') solves Lo·X=Y
%   exploiting known sparsity in Y, which is such that Y = Yf(obs,miss) and the
%   (i,j)-block of Yf is known to be zero for i > max(p,j+g). This function is
%   used to form Vhat (with g=0) and Laomh (with g=p) in the pure autoregressive
%   case. Lu, LSig and JPAT come from omega_ar_factor and give the ltl-factor of
%   Omega_o, ko(t) is the number of observed and km(t) the number of missing
%   values before time t. X = PROFILE_AR_TRISOLVE(...,'T') solves Lo'·X=Y.
%
%   [X,Xd] = PROFILE_AR_TRISOLVE(...,Lud,Lld,Yd) finds also the derivative of X.
%
%   NOTE: The routine is implemented using a blocked algorithm for speed

function [X,Xd] = profile_ar_trisolve(Lu,LSig,Y,JPAT,ko,km,p,g,cod,Lud,LSigd,Yd) 
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
        X(I,K) = omega_ar_trisolve(Lu, LSig, Y(I,K), JPAT, ko, cod);
      else
        [X(I,K),Xd(I,K,:)] = ...
          omega_ar_trisolve(Lu, LSig, Y(I,K), JPAT, ko,cod,Lud,LSigd,Yd(I,K,:));
      end
    end
  end
end
