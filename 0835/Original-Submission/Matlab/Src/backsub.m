
function x = backsub(A,b)
%
% assume A is an upper triangular matrix
% backsub performs backward substitution to solve Ax=b
%  syntax: >>x = backsub(A,b)
%
   n = size(A,2);
   x = zeros(n,1);
   
   B = diagproc(A);
   
   x(n) = b(n)/B(n,n);
   for k = (n-1):-1:1
       x(k) = (b(k) - B(k,k+1:n)*x(k+1:n))/B(k,k);
   end;
