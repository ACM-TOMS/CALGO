function A = sylves(f,g,k)
%
%  The k-th sylvester matrix of f(x), f'(x) = g(x)
%
   l = length(f);
   m = l+k-1;
   n = 2*k+1;
   A = zeros(m,n);
   for j = 1:k+1;
       A(j:j+l-2,j) = g';
   end;
   for j = 1:k
       A(j:j+l-1,j+k+1) = f';
   end;
