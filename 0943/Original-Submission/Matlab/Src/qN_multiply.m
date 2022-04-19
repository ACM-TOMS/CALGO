function [r] = qN_multiply(a,b,v,gamma_qn,stored_vecs)
%multiplies the L-BFGS matrix by v.
%initial Hessian =(1/gamma_qn)*I 

r = (1/gamma_qn)*v;

for i=1:stored_vecs
  r = r - (a(:,i)'*v)*a(:,i) + (b(:,i)'*v)*b(:,i);
end






