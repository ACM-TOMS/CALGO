
function z = hqrt(t,b)
%
%  transformation using the t from hessqr
%
   [m,n] = size(b);
   z = b;
   k = size(t,2);
   for j = 1:k
       T = [t(1,j),t(2,j); -t(2,j), t(1,j)];
       z(j:j+1,:)=T*z(j:j+1,:);
   end;
