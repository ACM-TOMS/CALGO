%THETA2MAT  Convert parameter vector to parameter matrices
%
%  [A,B,Sig] = THETA2MAT(THETA,p,q,r) converts a vector of the parameters of a
%  VARMA model to the parameter matrices A = [A1...Ap], B = [B1...Bq] and Sig.
%
%  [A,B,Sig,mu] = THETA2MAT(THETA,p,q,r) includes the mean mu among the
%  parameters.
%
%  [A,B,Sig] = THETA2MAT(THETA,p,q,r,J) allows for a linear change of variables,
%  so that [vec(A); vec(B); vech(Sig)] = J·THETA. If J has too few rows, an
%  identity matrix of the necessary size is appended to its bottom right: J :=
%  blkdiag(J, I), enabling, for instance, Sig to be excluded from the variable
%  change.
%
%  [A,B,Sig,mu] = THETA2MAT(THETA,p,q,r,J) sets [A(:); B(:); Sig(:)] to
%  J·THETA(1:k) and mu to THETA(k+1: end) where k is the column count of J.

function [A,B,Sig,mu] = theta2mat(theta,p,q,r,J)
  theta = theta(:);
  if nargin>=5
    m = size(J,2);
    theta = [J*theta(1:m); theta(m+1:end)];
  end
  i1 = p*r^2; i2 = i1 + q*r^2;
  A = reshape(theta(1:i1)   , r, r*p);
  B = reshape(theta(i1+1:i2), r, r*q);
  m = i2+1;
  for i=1:r
    m2 = m+r-i;
    Sig(i:r,i) = theta(m:m2);
    Sig(i,i:r) = theta(m:m2)';
    m = m2+1;
  end
  if nargout == 4, 
    if m<=length(theta), mu = theta(m:end); 
    else mu = []; end
  end
end
