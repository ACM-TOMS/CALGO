% CHNG_VAR  Postmultiply gradient with variable-change Jacobian
%
%    FD = CHNG_VAR(FD, J) changes derivatives FD from being w.r.t. A (and/or B
%    and Sig) to being w.r.t. an alternative set, theta, by postmultiplying the
%    leading columns of FD with the Jacobian of the variable change, J. See
%    help-text of var_ll and varma_llc.

function Fd = chng_var(Fd, J)
  if iscell(Fd)
    for i=1:length(Fd), Fd{i} = chng_var(Fd{i}, J); end
  else
    [mJ, nJ] = size(J);
    [a,b,np] = size(Fd);
    Fd1 = reshape(Fd(:,:,1:mJ),     a*b, mJ);
    Fd2 = reshape(Fd(:,:,mJ+1:end), a*b, np-mJ);
    Fd = reshape([Fd1*J Fd2], a, b, nJ + np-mJ);
  end
end
