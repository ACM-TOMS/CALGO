% ATA_DERIV  Calculate A'·A and its derivative
%
%   [C, Cd] = ATA_DERIV(A, Ad) calculates C = A'·A and the derivative of C with
%   respect to parameter theta(i) in Cd(:,:,i), i = 1,...,nPar. Ad(:,:,i) should
%   be the derivatives of A w.r.t. theta(i).
%
%   Use Cd = ATA_DERIV(A,Ad) when only the derivative is needed.

function [varargout] = ata_deriv(A, Ad)
  [m,n,nPar] = size(Ad);
  M = A'*reshape(Ad,m,n*nPar);
  Cd = reshape(M, n, n, nPar);
  Cd = Cd + permute(Cd,[2,1,3]);
  if nargout == 1, 
    varargout{1} = Cd;
  else
    varargout = {A'*A, Cd};
  end
end
