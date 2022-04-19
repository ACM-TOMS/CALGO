% LAMBDA_MULTIPLY  Multiply with Lambda
%
%  Y = LAMBDA_MULTIPLY(A, X, false(r,n)) returns Lam·X in Y, where Lam is an n·r
%  × n·r lower block-band matrix given by:
%
%               Lam =   I
%                           I
%                               I
%                                   I
%                      -Ap ... -A2 -A1  I
%                          -Ap     -A2 -A1  I
%                              ...         ..  ..
%                                  -Ap ...    -A1  I
%
%  and X is an (n·r)-element column vector. A should be [A1...Ap].
%
%  Y = LAMBDA_MULTIPLY(A, X, miss) where miss is r × n and true in missing
%  locations multiplies with Lam(obs,obs) where obs are the observed indices in
%  {1,2,...,r·n}. X should be an N-vector, where N is the number of observed
%  values. This call is used for multiplying with the matrix Lam_o.
%
%  [Y,Yd] = LAMBDA_MULTIPLY(A,X,miss,Xd) returns Lam·X in Y and the derivatives
%  of it in Yd, when Xd contains the derivatives of X w.r.t. all the parameters.
%  Yd is an N × 1 × nPar 3-dimensional array. 

function [Y, Yd] = lambda_multiply(A, X, miss, Xd)
  DIFF = nargout > 1;
  [r, n] = size(miss);
  MISS = any(miss(:));
  p = size(A,2)/r;
  Ac = makecell(A);
  if MISS
    obs = ~miss;
    Xf = zeros(r*n,1);
    Xf(obs) = X;
    X = reshape(Xf, r, n);
    [ko, ro, km] = find_missing_info(miss);
    N = ko(n+1);
  else
    X = reshape(X, r, n);
  end
  Y = X;
  if DIFF
    nPar = size(Xd,3);
    if MISS
      Xfd = zeros(n*r, nPar);
      Xfd(obs,:) = Xd;
    else
      Xfd = Xd;
    end
    I = any(Xfd,1); % nonzeros in nPar direction
    nI = sum(I);
    Xfd = reshape(Xfd, r, n, nPar);
    Zd = permute(Xfd(:,:,I),[1,3,2]);
    Xd = reshape(Zd, r*nI, n);
    L = p*nPar+1 : n*nPar;
  end
  L = p+1:n;
  for k = 1:p
    K = p+1-k:n-k;
    Y(:, L) = Y(:, L) - Ac{k}*X(:, K);
    if (DIFF)
      YdIL = reshape(Zd(:,:,L), r, nI*(n-p));
      YdIL = YdIL - Ac{k}*reshape(Xd(:,K), r, nI*(n-p));
      Zd(:,:,L) = reshape(YdIL, r, nI, n-p);
    end
  end
  Y = Y(:);
  if MISS, Y = Y(obs); end
  if (DIFF)
    Yd = zeros(r, n, nPar);
    Yd(:,:,I) = permute(reshape(Zd, r, nI, n), [1,3,2]);
    for k=1:p
      j0 = r^2*(k-1);
      for j = 1:r
        J = j0 + (j:r:r^2); lJ = length(J);
        Yd(j, p+1:n, J) = Yd(j, p+1:n, J) - reshape(X(:,p-k+1:n-k)', 1, n-p,lJ);
      end
    end
    Yd = reshape(Yd, n*r, 1, nPar);
    if MISS, Yd = Yd(obs, 1, :); end
  end
end
