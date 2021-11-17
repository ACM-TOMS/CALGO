function [B] = qN_B(a,b,gamma_qn,stored_vecs)
%computes B, the LBFGS matrix
%initial Hessian =(1/gamma_qn)*I 
%only used by ms.m

B = (1/gamma_qn)*eye(size(a,1));

for i=1:stored_vecs
  B = B - a(:,i)*a(:,i)' + b(:,i)*b(:,i)';
end






