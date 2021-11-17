%MDS_SET_PARMAT  Set matrix and derivative to a parameter matrix and derivative
%
%  FS = MDS_SET_PARMAT(nparmat,X,kX) sets FS to a structure appropriate for
%  mds_add_prod. FS.mat is set to X which should be the kX-th paramter matrix
%  (out of the sequence A1...Ap B1...Bq Sig), and nparmat=p+q+1. FS.der is set
%  to the derivative of X with respect to all the parameter matrices
%  FS.spcod(kX) is set to 'e' and FS.spcod(j) is set to 'z' for j~=kX.

function FS = mds_set_parmat(nparmat,X,kX)
  r = size(X,1);
  FS = mds_set_zero(nparmat,r,r);
  FS.mat = X;
  for c=1:r
    for l=1:r
      FS.der{kX}{l,c}(l,c) = 1;
    end
  end
  FS.spcod(kX) = 'e';
end
