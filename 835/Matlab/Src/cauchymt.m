
function A = cauchymt(f,k)
   %
   % program to generate a k-column Cauchy matrix of polynomial f
   %
   n = length(f);
   m = n + k - 1;
   A = zeros(m,k);
   
   for j = 1:k
       A(j:j+n-1,j) = f.';
   end;
