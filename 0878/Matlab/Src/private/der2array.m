%DER2ARRAY  Convert derivative from cell format to array format
%
%  [F,Fd] = DER2ARRAY(FS) extracts the cell arrays F and Fd from the struct
%  array FS. F{i} is r×r and Fd{i} is r×r×NPAR with Fd{i}(:,:,m) containing the
%  derivative of F{i} with respect to the m-th parameter; NPAR = r^2·(p+q) +
%  r·(r+1)/2.
%
%  The parameters are ordered by parameter matrix, A1,..., Ap, B1,..., Bq, Sig
%  and within each parameter matrix by columns, with the upper triangle of Sig
%  removed: A1(1,1), A1(2,1),..., Bq(r,r), Sig(1,1), Sig(r,1), Sig(2,2),...,
%  Sig(r,r).
%
%  Apart from simple copying, DER2ARRAY adjusts for the fact that Sig is
%  symmetric. Derivatives w.r.t. mu (if present) are set to zero. FS is a
%  structure array, the i-th element containing F{i} and its derivatives in a
%  format as described in mds_add_prod.
%
%  Fd = der2array(FS) returns only the derivative.

function [F, Fd] = der2array(FS)
  n = length(FS);
  r = size(FS(1).mat, 1);
  nPar = (length(FS(1).der)-1)*r^2 + r*(r+1)/2;
  F = cell(1,n);
  if nargout>1
    Fd = cell(1,n);
    for i = 1:n, [F{i},Fd{i}] = dercopy(FS(i), nPar); end
  else
    for i = 1:n, F{i} = dercopy(FS(i), nPar); end
  end
end

function [F, Fd] = dercopy(FS, nPar);
  r = size(FS.mat,1);
  nparmat = length(FS.der);
  Fd = zeros(r,r,nPar);
  m = 0;
  for k=1:nparmat-1 % Copy derivatives w.r.t. A and B matrices
    for c=1:r
      for l=1:r
        m = m+1;
        Fd(:,:,m) = FS.der{k}{l,c};
      end
    end
  end
  for c=1:r % Copy derivative w.r.t. Sigma, and adjust for symmetry
    m = m+1;
    Fd(:,:,m) = FS.der{nparmat}{c,c};
    for l=c+1:r
      m = m+1;
      Fd(:,:,m) = FS.der{nparmat}{l,c} + FS.der{nparmat}{c,l};
    end
  end
  ascertain(m==nPar || m==nPar-r);
  if nargout>1
    F = FS.mat;
  else
    F = Fd;
  end
end
