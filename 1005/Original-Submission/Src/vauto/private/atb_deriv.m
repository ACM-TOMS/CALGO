%ATB_DERIV  Calculate A'·B and its derivative
%
%  [C, Cd] = ATB_DERIV(A, Ad, B, Bd) calculates C = A'·B and the derivative of C
%  with respect to parameter theta(i) in Cd(:,:,i), i = 1,...,nPar. Ad(:,:,i)
%  and Bd(:,:,i) should be the derivatives of A and B w.r.t. theta(i).
%
%  Use Cd = ATB_DERIV(A, Ad, B, Bd) when only the derivative is needed.
%
%  METHOD: Cd(:,:,i) = Ad(:,:,i)'*B + A'*Bd(:,:,i).

function [varargout] = atb_deriv(A, Ad, B, Bd)
  [k,m] = size(A);
  n = size(B,2);
  nPar = size(Ad,3);
  BA = reshape(B'*reshape(Ad,k,m*nPar), n, m, nPar);
  BA = permute(BA, [2,1,3]);
  Cd = reshape(A'*reshape(Bd,k,n*nPar), m, n, nPar) + BA;
  if nargout == 1, 
    varargout{1} = Cd;
  else
    [varargout{1:2}] = deal(A'*B, Cd);
  end
end
