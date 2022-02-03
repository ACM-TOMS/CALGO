function [a,b]=unrolling(b0)

global S Y RHO;

% "Unrolling the BFGS formula" 
%
% This code is an implementation of Procedure 7.6 found in 
% J. Nocedal, S.J. Wright, Numerical Optimization, Springer-Verlag, New
% York, second edition, 2006.  

% Inputs: S, Y contain the L-BFGS updates 
% where B_0 = b0*I is the initial L-BFGS matrix

% The most recently computed s_i and y_i vectors are stored in the last
% columns of S and Y, respectively.

[m,n] = size(S);
a = zeros(m,n); 
b = zeros(m,n); 


for i = 1:n
  b(:,i) = Y(:,i)*sqrt(RHO(i)); %Y(:,i)/(sqrt(Y(:,i)'*S(:,i)));  
  a(:,i) = b0*S(:,i);
  for j = 1: (i-1)
    a(:,i) = a(:,i) + b(:,j)'*S(:,i)*b(:,j) - a(:,j)'*S(:,i)*a(:,j);
  end
  a(:,i) = a(:,i)/(sqrt(S(:,i)'*a(:,i)));
end



