% FORWARD_SUB_DERIV  Derivative of forward substitution
%
%  [X, Xd] = FORWARD_SUB_DERIV(L, Ld, Y, Yd) solves a lower triangular system
%  L·X=Y and returns the derivative of X in Xd.
%
%  Xd = FORWARD_SUB_DERIV(L, Ld, X, Yd) finds the derivative if X is given.
%
%  The parameters are:
%    L   n×n        A lower triangular matrix
%    Ld  n×n×nPar   Ld (:,:,i) contains derivatives of L with respect to a
%                   parameter theta(i), i=1...,nPar
%    X   n×m        Solution to L·X = Y
%    Yd  n×m×nPar   Derivatives of Y. May also be n×nPar when m=1
%    Xd  n×m×nPar   Returns derivatives of X. Is n×nPar if Yd is n×nPar.
%
%  METHOD:
%    From L·X = Y it follows that Ld(:,:,i)·X + L·Xd(:,:,i) = Yd, and thus Xd
%    may be obtained with forward substitution. 
%
function [varargout] = forward_sub_deriv(L, Ld, varargin)
  LT.LT = true;
  if nargout == 2
    [Y, Yd] = deal(varargin{:});
    X = linsolve(L, Y, LT);
  else
    [X, Yd] = deal(varargin{:});
  end
  [n, m] = size(X);
  nPar = size(Ld, 3);
  LdX= reshape(reshape(permute(Ld, [1,3,2]), n*nPar, n)*X, n,nPar, m);
  RHS = reshape(Yd - permute(LdX, [1,3,2]), [n,m*nPar]);
  Xd = reshape(linsolve(L,RHS,LT), [n,m,nPar]);
  if nargout == 2, varargout = {X, Xd}; else varargout = {Xd}; end
end
