% VYW_DERIV_RHS  Rhs of equations for derivative of vector-Yule-Walker solution
%
%    RHS = VYW_DERIV_RHS(A,GGd,S) calculates the right-hand-side of equations
%    that determine the derivative of the Sj matrices calculated by VYW. A is a
%    cell array of the AR parmeter matrices, GGd is an mds struct-array of Gi
%    derivatives (calculated with find_CGW_deriv), and S should be previously
%    calculated with vyw_solve. The i-th right-hand-side matrix is returned in
%    the der components of the r×r mds structure RHS{i}, i=1,...,p+1.

function RHS = vyw_deriv_rhs(A, GGd, S)
  r = size(GGd(1).mat,1);
  p = size(A,2)/r;
  q = length(GGd) - 1;
  A = makecell(A);
  RHS = cell(1,p+1);
  for j=0:p
    if j<=q, R = GGd(j+1); end
    if j> q, R = mds_set_zero(p+q+1,r,r); end
    for i = 1:j
      R = mds_add_prod(R, A{i}, S{j-i+1}, i);
    end
    for i = j+1:p
      R = mds_add_prod(R, A{i}, S{i-j+1}', i);
    end
    RHS(j+1:j+1) = der2array(R);
  end
end
