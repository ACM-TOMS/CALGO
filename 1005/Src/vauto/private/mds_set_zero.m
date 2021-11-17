%MDS_SET_ZERO  Set matrix and its derivative to zero
%
%  FS = MDS_SET_ZERO(p+q+1,r,n) sets FS to a structure appropriate for
%  mds_add_prod containing an r×n matrix and its derivative both equal to 0; p,
%  q and r are the dimensions of the VARMA model being applied.
%
function FS = mds_set_zero(npar,r,n)
  FS.mat = zeros(r,n);
  FS.der = cell(1,npar);
  der_k = mat2cell(zeros(r^2,n*r),repmat(r,r,1),repmat(n,r,1));
  for k = 1:npar
    FS.der{k} = der_k;
  end
  FS.spcod = repmat('z',1,npar);
end
