%CHOL_DERIV  Derivative of Cholesky factors
%
%  Ld = CHOL_DERIV(L, Ad) is used to differentiate Cholesky factorization. If A
%  = L*L' is the Cholesky factorization of an n×n matrix A and Ad(:,:,i)
%  contains derivatives of A with respect to parameter theta(i) for i=1,...,nPar
%  then Ld(:,:,i) will be the derivative of L with respect to theta(i). Only the
%  lower triangle of Ad(:,:,i) is used.
%
%  [L, Ld] = CHOL_DERIV(A, Ad) finds the Cholesky factor of A in L and its
%  derivative in Ld.

function [varargout] = chol_deriv(varargin)
  LT.LT = true;
  if nargout == 2
    [A, Ad] = deal(varargin{:});
    L = cholf(A);
  else
    [L, Ad] = deal(varargin{:});
  end
  n = length(L);
  nPar = size(Ad,3);
  Ld = zeros(n, nPar*n);
  for k=1:n
    g = nPar*(k-1);
    K = 1:k-1;
    u = L(k,K);
    X = zeros(k-1,nPar);
    blk = 15; % use blocking to exploit fast multiplication for larger matrices
    j = 1;
    for i=1:blk:k-1
      m = min(blk, k-i);
      ie = i+m-1;
      je = j+m*nPar-1;
      X(i:ie,:) = reshape(u(1:ie)*Ld(1:ie,j:je), nPar, m)';
      j = je+1;
    end
    % X(:,i) now contains Ldi(K,K)*u'
    X1 = reshape(Ad(k,K,:),k-1,nPar);
    X = linsolve(L(K,K), X1-X, LT);
    Ld(K, g+1:g+nPar) = X;
    sd = reshape(Ad(k,k,:),1,nPar)/2;
    Ld(k, g+1:g+nPar) = (sd - u*X)/L(k,k);
  end
  Ld = permute(reshape(Ld, n, nPar, n), [3,1,2]);
  if nargout == 2, varargout = {L, Ld}; else varargout = {Ld}; end
end
