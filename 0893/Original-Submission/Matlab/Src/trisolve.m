function x = trisolve(a,b,c,d)
% trisolve:  Solution to nonsymmetric tridiagonal linear system
%
% USAGE:  x = trisolve(a,b,c,d);
%
% This function computes the solution to the equations
% 
%    b(1)*x(1) + c(1)*x(2) = d(1)
%    a(i-1)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i), i = 2:n-1
%    a(n-1)*x(n-1) + b(n)*x(n) = d(n)
%
% No pivoting is used.
%
% On input:
%
%       A,B,C = Vectors with lengths n-1, n, and n-1, 
%               respectively, containing the subdiagonal,
%               diagonal, and superdiagonal of the matrix.
% 
%       D = Array dimensioned n by k containing one or more
%           right hand side vectors.
%
% On output:
%
%       X = Array of size [n k] containing the solution 
%           vectors.
%
% Modules required by TRISOLVE:  None
%
%***********************************************************

n = length(b);
x = d;

for i = 1:n-1                        % Forward elimination
   mu = -a(i)/b(i);
   b(i+1) = b(i+1) + mu*c(i);
   x(i+1,:) = x(i+1,:) + mu*x(i,:);
end

x(n,:) = x(n,:)/b(n);                % Back substitution
for i = n-1:-1:1
   x(i,:) = (x(i,:)-c(i)*x(i+1,:))/b(i);
end
return;

end  % trisolve
