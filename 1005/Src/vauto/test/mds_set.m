%MDS_SET_PARMAT  Set matrix and derivative to given values
%
%  FS = MDS_SET(F,Fd) sets FS to a structure appropriate for mds_add_prod.
%  FS.mat is set to F, FS.der{k} is set to Fd(:,:,r^2·(k-1)+1:r^2·k) and
%  FS.spcod(k) is set to 'f' for all k. The operation is the inverse of
%  der2array.

function FS = mds_set(F,Fd)
  r = size(F,1);
  nPar = size(Fd,3);
  nparmat = (nPar + r*(r-1)/2)/(r*r);
  FS = mds_set_zero(nparmat,r,r);
  FS.mat = F;
  FS.spcod(nparmat) = 'f';
  m = 0;
  for k = 1:nparmat-1
    for c=1:r
      for l=1:r
        m = m+1;
        FS.der{k}{l,c} = Fd(:,:,m);
      end
    end
    FS.spcod(k) = 'f';
  end
  for c=1:r % Copy derivative w.r.t. Sigma, and adjust for symmetry
    m = m+1;
    FS.der{nparmat}{c,c} = Fd(:,:,m);
    for l=c+1:r
      m = m+1;
      FS.der{nparmat}{l,c} = Fd(:,:,m)/2;
      FS.der{nparmat}{c,l} = Fd(:,:,m)/2;
    end
  end
end
