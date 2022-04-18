function [U,R] = sylup(f,g,k,U,R)
   % 
   % updating the Sylvester discriminant matrix and QR decomposition
   %
   n = length(f);
   if k == 1 
       %
       % construct S_1
       %
       S = zeros(n,3);  S(1:n-1,1) = g.'; S(:,2) = f.';  S(2:n,3) = g.';
       %
       % triangularize it
       %
       U = zeros(n,3); R = S;
       for j = 1:3
           if R(j,j) ~= 0
               s = norm(R(j:n,j))*R(j,j)/abs(R(j,j));
           else
               s = norm(R(j:n,j));
           end;
           u = R(j:n,j); u(1) = u(1) + s; u = u/norm(u);
           U(j:n,j) = u;
       
           R(j:n,j:3) = R(j:n,j:3) - 2*u*(u'*R(j:n,j:3));
       end;
   else
       %
       % extend S_j
       %
       R = [R; zeros(1,2*k-1)]; U = [U; zeros(1,2*k-1)];
       R = [R, [zeros(k-1,1);f.'],[zeros(k,1);g.']];
       U = [U, zeros(n+k-1,2)];
       %
       % update previous transformation
       %
       jj = [2*k, 2*k+1]; 
       for j = 1:2*k-1
           m = n+max(0,floor((j-2)/2));
           u = U(j:m,j);
           R(j:m,jj) = R(j:m,jj) - 2*u*(u'*R(j:m,jj));
       end;
       %
       % additional transformations
       %
       m = n+k-1;
       for j = 2*k:2*k+1
           if R(j,j) ~= 0
               s = norm(R(j:m,j))*R(j,j)/abs(R(j,j));
           else
               s = norm(R(j:m,j));
           end;
           u = R(j:m,j); u(1) = u(1) + s; 
           t = norm(u); if t ~= 0, u = u/t; end;
           U(j:m,j) = u;
       
           R(j:m,j:2*k+1) = R(j:m,j:2*k+1) - 2*u*(u'*R(j:m,j:2*k+1));
       end;
           
   end;
