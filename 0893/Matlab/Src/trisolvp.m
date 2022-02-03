function x = trisolvp(a,b,c,d,l,u)
% trisolvp:  Solution to nonsymmetric almost tridiagonal linear system
%
% USAGE:  x = trisolvp(a,b,c,d,l,u);
%
% This function computes the solution to the equations
% 
%    b(1)*x(1) + c(1)*x(2) + u*x(n) = d(1)
%    a(i-1)*x(i-1) + b(i)*x(i) + c(i)*x(i+1) = d(i), i = 2:n-1
%    l*x(1) + a(n-1)*x(n-1) + b(n)*x(n) = d(n)
%
% The matrix is tridiagonal with nonzeros in the upper 
% right and lower left corners.  This arises from periodic
% end conditions.  No pivoting is used.
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
%       L,U = Lower left and upper right corner elements of
%             the matrix.
%
% On output:
%
%       X = Array of size [n k] containing the solution 
%           vectors.
%
% Modules required by TRISOLVP:  None
%
%***********************************************************

n = length(b);
u = [u; zeros(n-2,1)];  % Extend length(u) to allow fill-in.
x = d;

for i = 1:n-2                          % Forward elimination
   mu = -a(i)/b(i);
   b(i+1) = b(i+1) + mu*c(i);
   u(i+1) = mu*u(i);
   x(i+1,:) = x(i+1,:) + mu*x(i,:);
   mu = -l/b(i);
   l = mu*c(i);
   b(n) = b(n) + mu*u(i);
   x(n,:) = x(n,:) + mu*x(i,:);
end
u(n-1) = u(n-1) + c(n-1);
l = l + a(n-1);
mu = -l/b(n-1);
b(n) = b(n) + mu*u(n-1);

x(n,:) = (x(n,:) + mu*x(n-1,:))/b(n);    % Back substitution
x(n-1,:) = (x(n-1,:)-u(n-1)*x(n,:))/b(n-1);
for i = n-2:-1:1
   x(i,:) = (x(i,:)-c(i)*x(i+1,:)-u(i)*x(n,:))/b(i);
end
return;

end  % trisolvp
