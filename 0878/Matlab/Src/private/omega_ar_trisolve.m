% OMEGA_AR_TRISOLVE Triangular solve with Cholesky factor of Omega, pure AR case
%
%   X = OMEGA_AR_TRISOLVE(Lu, LSig, Y, JPAT, ko, 'NoT') solves L·X = Y where L
%   is the Cholesky factor of Omega(obs,obs) in a pure auto-regressive case
%   (q=0) stored in Lu and LSig, which are returned by omega_ar_factor. ko(t) is
%   the number of observed values before time t.
%
%   X = OMEGA_AR_TRISOLVE(Lu, LSig, Y, JPAT, ko, 'T') solves L'·X = Y.
%
%   [X, Xd] = OMEGA_AR_TRISOLVE(Lu, LSig, Y, JPAT, ko, code, Lud, LSigd, Yd)
%   finds also the deriative of X. Lud, LSigd and Yd are the derivatives of Lu,
%   LSig and Y, and code begins with 'N' or 'T'
%
%   As in omega_back_sub, Y may have fewer than N rows (N=number of observed
%   values), and then the rest of the rows are taken to be zero. It must end
%   on a block boundary. 

function [Y,Yd] = omega_ar_trisolve(Lu, LSig, Y, jpat, ko, code, Lud, LSigd, Yd)
  % Note: Xd and Xd are stored in Y and Yd
  OPTS.LT = true;
  if code(1)=='T', OPTS.TRANSA = true; end
  n = length(ko)-1;
  p = n - length(jpat);
  mY = size(Y,1);
  tmax = find(ko == mY, 1, 'first') - 1; ascertain(~isempty(tmax));
  e = size(Lu,1);
  [mY,nY] = size(Y);
  Y(1:e, :) = linsolve(Lu, Y(1:e, :), OPTS);
  for t=p+1:tmax
    K = ko(t)+1 : min(mY, ko(t+1));
    Y(K,:) = linsolve(LSig{jpat(t-p)}, Y(K,:), OPTS);
  end  
  if nargout > 1,
    nPar = size(Yd, 3);
    Lud = reshape(permute(Lud, [1,3,2]), e*nPar, e);
    Yd(1:e,:,:) = Yd(1:e,:,:)-permute(reshape(Lud*Y(1:e,:),e,nPar,nY), [1,3,2]);
    if isempty(Yd), Xd = Yd; return; end
    for t=p+1 : tmax
      K = ko(t)+1 : min(mY, ko(t+1));
      ro = ko(t+1) - ko(t);
      Lj = reshape(permute(LSigd{jpat(t-p)}, [1,3,2]), ro*nPar, ro);
      if code(1)=='T'
        Yd(K,:,:) = Yd(K,:,:) - permute(reshape(Lj'*Y(K,:),ro,nPar,nY),[1,3,2]);
      else
        Yd(K,:,:) = Yd(K,:,:) - permute(reshape(Lj*Y(K,:),ro,nPar,nY),[1,3,2]);
      end
    end
    Yd(1:e, :) = linsolve(Lu, Yd(1:e, :), OPTS);
    for t=p+1:tmax
      K = ko(t)+1 : min(mY,ko(t+1));
      Yd(K,:) = linsolve(LSig{jpat(t-p)}, Yd(K,:), OPTS);
    end   
    Yd = reshape(Yd, mY, nY, nPar);
  end
end
