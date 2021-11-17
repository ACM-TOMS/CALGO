%FIND_CGW  Determine C, G and W matrices for VARMA log-likelihood
%
%  [C,G,W] = FIND_CGW(A,B,Sig) calculates the matrices Ci, Gi and Wi for
%  VARMA_LLC (see comments there). On entry A={A1 A2...Ap} and B={B1 B2...Bq}
%  are cell arrays of length p and q with the r×r matrices Ai and Bi and Sig is
%  also r×r. On exit C, G and W are cell arrays of r×r matrices:
%                C = {C0 C1...Cq} with Cj = cov(x(t),eps(t-j)),
%                G = {G0 G1...Gq} with Gj = cov(y(t),x(t-j)), and 
%                W = {W0 W1...Wq} with Wj = cov(y(t),y(t-j)).
%  These matrices are given by the formulae:
%                Cj = A1·C(j-1) + ... + A(j-1)·C1 + Aj·C0 + Bj·Sig 
%                Gj = Bj·C0' + ... + Bq·C(q-j)'
%                Wj = Bj·Sig·B0' + ... + Bq·Sig·B(q-j)'
%  where C0 = Sig and B0 = -I.

function [C,G,W] = find_CGW(A,B,Sig)
  [p,q,r] = get_dimensions(A,B,Sig);
  I = eye(r);
  A = makecell(A);
  B = [{I} makecell(B)];
  C = cell(1,q+1);
  G = cell(1,q+1);
  W = cell(1,q+1);
  BSig = cell(1,q+1);
  C{1} = Sig;
  for j=0:q
    BSig{j+1} = B{j+1}*Sig;
    C{j+1} = BSig{j+1};
    for i=1:min(j,p)
      C{j+1} = C{j+1} + A{i}*C{j-i+1};
    end
  end
  for j = 0:q
    G{j+1} = BSig{j+1};
    W{j+1} = BSig{j+1};
    for i = j+1:q
      G{j+1} = G{j+1} + B{i+1}*C{i-j+1}';
      W{j+1} = W{j+1} + BSig{i+1}*B{i-j+1}';
    end
  end
end
