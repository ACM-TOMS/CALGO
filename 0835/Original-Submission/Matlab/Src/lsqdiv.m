
function u = lsqdiv( f, v, w)
%
%  polynomial division u = f/v using least squares method
%  generate a k-column Cauchy matrix of polynomial f
%
   m = length(f);     k = length(v);     n = m - k + 1;
   A = zeros(m,n);
   
   if nargin == 2   % construct weights, if not provided
       w = ones(m,1);
       for j = 1:m,  if abs(f(j)) > 1, w(j) = 1/abs(f(j)); end;  end
   end;
   %
   % construct the augmented Cauchy matrix
   %
   kk = k-1;   
   for j = 1:n,      A(j:j+kk, j) = v.';       end;
   A = [A,f.'];
   for j = 1:m,     A(j,:) = A(j,:)*w(j);    end;  % scale
   %
   % QR decomposition of the Cauchy matrix
   %
   U = [];
   for j = 1:n
       x = A(j:j+kk, j);
       s = norm(x);  if x(1) ~= 0, s = s*x(1)/abs(x(1));  end;
       u = x;              u(1) = u(1) + s;      u = u/norm(u);
       A(j,j) = -s;        A(j+1:j+kk, j)  =  zeros(kk,1);
       ll = [j+1:min(j+kk,n), n+1];
       for l = ll
           s = 2*u'*A(j:j+kk, l);
           A(j:j+kk, l) = A(j:j+kk, l) - s*u;
       end;
       U = [U,u];
   end;
   %
   % banded backward substitution
   %
   u = zeros(n,1); b = A(1:n,n+1);
   u(n) = b(n)/A(n,n);
   %
   % triangular part
   %
   for j = (n-1):-1:max(1,n-kk)
       u(j) = (b(j) - A(j,j+1:n)*u(j+1:n))/A(j,j);
   end;
   %
   % banded part
   %
   if n-k >=1 
       for j = n-k:-1:1
           u(j) = (b(j)-A(j,j+1:j+k)*u(j+1:j+k))/A(j,j);
       end;
   end;
