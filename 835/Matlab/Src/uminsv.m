
function [s,x] = uminsv(R,tol,nit)
%
% Find mininum singular value of a triangular matrix R
%
   if nargin == 2, nit = 3; end;
   [m,n] = size(R);      % get the dimensions of R
   scale = norm(R,inf);  % get the magnitude of rows as scaler
   a = 2*rand(1,n)-1;        % random initial vector
   a = scale*a/norm(a);      % set the first row
   b = [scale;zeros(m,1)];
     
   [T,trans] = hessqr([a;R]);  % Hessenberg QR decomp. of stacked matrix
   z = hqrt(trans,b);          % same Q on b
   
   x = backsub(T(1:n,1:n),z(1:n));  x = x/norm(x); % geting the new vector
   
   r = [scale*x';R]*x-b;
   ss = [norm(r(2:m+1))]; cr = [];
   
   for k = 1:nit
       
       [T,trans] = hessqr([2*scale*x';R]); 
       z = hqrt(trans,r); 
       u = backsub(T(1:n,1:n),z(1:n));
       y = x - u;
       y = y / norm(y);
       r = [scale*y';R]*y-b; 
       s = norm(r(2:m+1));
       ss = [ss,s];
       cr = [cr,norm(x-y)]; %disp(cr(k));
      
       if k == 1
           if cr(1) < tol, break; end;
       else
           if cr(k) < cr(k-1)
               if cr(k)^2/(cr(k-1)-cr(k))<tol 
                   break
               end;
           end;
       end;
       x = y;
       
   end;
