
function x = scalsq(A,b,w)
%
%  Solving least squares problem with iterative refinement
%      
   [m,n] = size(A);
   
   if nargin == 2
       w = ones(m,1);
       for j = 1:m
           if abs(b(j)) > 1, w(j) = 1/abs(b(j)); end;
       end
   end;
   
   for j = 1:m
       A(j,:) = A(j,:)*w(j);
   end;
   b = b.*w;
   
   [Q,R] = qr(A);

   d = Q'*b;  S = R(1:n,1:n);
   x = backsub(S,d(1:n));
   
   % one step refinement
   bb = [b;zeros(n,1)]; B = [eye(m),A;A',zeros(n,n)];
   r = b - A*x;
   
   %for j = 1:3
   rr = bb - B*[r;x]; %disp(norm(rr));
   
   s = Q'*rr(1:m);
   c = forsub(S',rr(m+1:m+n));
   c2 = backsub(S,s(1:n)-c);
   c1 = Q*[c;s(n+1:m)];
   %r = r + c1;
   x = x + c2;
   %end;
