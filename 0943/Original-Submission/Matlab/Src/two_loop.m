function [r] = two_loop(vec,gamma_qn,stored_vecs)
%multiplies "vec" with the inverse of the L-BFGS matrix
%assumes initial L-BFGS inverse is  gamma_qn*I, i.e., B_0=(1/gamma_qn)I 
%an implementation based on the two-loop recursion published in
%NOCEDAL, J. 1980. Updating quasi-Newton matrices with limited storage. 
%Math. Comput. 35, 773--782.

global S Y RHO;
q = vec;
for i = stored_vecs:-1:1
  alfa(i) = RHO(i)*S(:,i)'*q;
  q = q - alfa(i)*Y(:,i);
end
r = gamma_qn*q;
for i = 1:stored_vecs
  beta = RHO(i)*Y(:,i)'*r;
  r = r + (alfa(i)-beta)*S(:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


