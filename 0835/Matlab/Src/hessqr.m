
function [B,t] = hessqr(A)
%
%  assume A is a Hessenberg matrix
%  output  B -- upper triangular
%          c -- rotation used
%
   [m,n] = size(A);
   if m < n, return; end;
   
   B = A;
   
   for j = 1:n
       if j < m
           d = sqrt(B(j,j)^2+B(j+1,j)^2);
           c = B(j,j)/d; s = B(j+1,j)/d;
           T = [c,s;-s,c];
           B(j:j+1,j:n) = T*B(j:j+1,j:n);
           t(1,j) = c; t(2,j) = s;
       end;
   end;
