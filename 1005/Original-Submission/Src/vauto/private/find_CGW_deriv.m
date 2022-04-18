%FIND_CGW_DERIV  Determine G, W and their derivatives
%
%  [CCd, GGd, WWd] = find_CGW_deriv(A, B, Sig) returns three struct arrays, CCd
%  = [C0...Cq], GGd = [G0...Gq] and WWd = [W0...Wq] where each element is a
%  structure of a matrix and its derivatives w.r.t. all the parameters as
%  described in mds_add_prod. Thus CCd(t).mat contains Ct and CCd(t).der{k}{l,c}
%  contains the derivatives of Ct w.r.t. the (l,k) element of the k-th parameter
%  matrix (of A1..Ap B1..Bq Sig).
%
function [C, G, W] = find_CGW_deriv(A, B, Sig)
  [p, q, r] = get_dimensions(A, B, Sig);
  A = makecell(A);
  B = makecell(B);
  np = p+q+1;
  for j=1:q
    BSig(j) = mds_set_zero(np, r, r);
    BSig(j) = mds_add_prod(BSig(j), B{j}, Sig, p+j, p+q+1);
    C(j) = BSig(j);
    for i=1:min(j-1, p)
      C(j) = mds_add_prod(C(j), A{i}, C(j-i), i);
    end
    if j <= p, C(j) = mds_add_prod(C(j), A{j}, Sig, j, p+q+1); end
  end
  C0 = mds_set_parmat(np, Sig, np);
  G0 = mds_set_parmat(np, Sig, np);
  W0 = mds_set_parmat(np, Sig, np);
  for i = 1:q
    G0 = mds_add_prod(G0, B{i}, C(i), 'T', p+i); 
    W0 = mds_add_prod(W0, B{i}, BSig(i), 'T', p+i);
  end
  for j = 1:q
    G(j) = BSig(j);
    W(j) = BSig(j);
    for i = j+1:q
      G(j) = mds_add_prod(G(j), B{i}, C(i-j), 'T', p+i);
      W(j) = mds_add_prod(W(j), B{i}, BSig(i-j), 'T', p+i);
    end
  end
  if q==0, C=[]; G=[]; W=[]; end
  C = [C0, C];
  G = [G0, G];
  W = [W0, W];
end
