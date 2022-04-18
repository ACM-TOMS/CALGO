%AAT_DERIV  Calculate A'·A and its derivative
%
%  [C, Cd] = AAT_DERIV(A, Ad) calculates C = A·A' and the derivative of C with
%  respect to parameter theta(i) in Cd(:,:,i), i = 1,...,nPar. Ad(:,:,i) should
%  be the derivatives of A w.r.t. theta(i).
%
%  Use Cd = AAT_DERIV(...) when only the derivative is needed.
%
function [varargout] = aat_deriv(A, Ad)
  [m,n] = size(A);
  nPar = size(Ad,3);
  Cd = zeros(m,m,nPar);
  for i = 1:nPar
    M = A*Ad(:,:,i)'; Cd(:,:,i) = M + M';
  end
  if nargout == 1, 
    varargout{1} = Cd;
  else
    [varargout{1:2}] = deal(A*A', Cd);
  end
end
