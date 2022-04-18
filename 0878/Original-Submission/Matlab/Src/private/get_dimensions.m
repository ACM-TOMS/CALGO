%DIMENSIONS  Dimensions of VARMA model
%
%  [p,q,r] = GET_DIMENSIONS(A, B, Sig) returns the dimensions of the VARMA
%  model with parameters A = {A1...Ap}, B = {B1...Bp} and Sig (Ai are the AR
%  coefficients, Bi are the MA coefficients and Sig is the covariance matrix of
%  the shock series.
%
%  [p,q,r,n] = GET_DIMENSIONS(A, B, Sig, X) returns also the length of the
%  time series (the number of rows in X) in n.

function [p, q, r, n] = get_dimensions(A, B, Sig, X)
  r = size(Sig,1);
  if r== 0, error('Sigma must not be empty'), end
  if iscell(A), p = length(A); else p = size(A,2)/r; end
  if iscell(B), q = length(B); else q = size(B,2)/r; end
  if nargin>3
    ascertain(size(X,1)==r);
    n = size(X,2); 
  end
end
