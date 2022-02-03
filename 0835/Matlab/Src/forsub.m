
function x = forsub(A,b)
%
% assume A is an lower triangular matrix
% forsub performs forward substitution to solve Ax=b
%
   n = size(A,1);
   x = zeros(n,1);
   
   B = diagproc(A);
   
   x(1) = b(1)/B(1,1);
   for k = 2:n
       x(k) = (b(k) - B(k,1:k-1)*x(1:k-1))/B(k,k);
   end;
