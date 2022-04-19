%MAT2THETA  Convert parameter matrices to a vector
%
%  THETA = MAT2THETA(A,B,Sig) changes the parameter matrices of a VARMA model, A
%  = {A1...Ap}, B = {B1...Bq} and Sig to a vector of all the parameters. It is a
%  utility used by functions that compare analytical derivatives with numerical
%  ones. A matrix input A = [A1 ... Ap] and B = [B1 ... Bq] is also permitted.
%
%  THETA = MAT2THETA(A,B,Sig,MU) adds MU at the end THETA.

function theta = mat2theta(A,B,Sig,mu)
  r = size(Sig,1);
  if iscell(A), A = cell2mat(A); end
  if iscell(B), B = cell2mat(B); end
  s = [];
  for i=1:r, s = [s; Sig(i:r,i)]; end
  theta = [A(:); B(:); s];
  if nargin==4, theta = [theta; mu];
end
