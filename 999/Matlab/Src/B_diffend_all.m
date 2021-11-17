function K = B_diffend_all(P, r, el)

% Full differentiation of all B-splines at one end point
% up to a given order
%
% INPUT
%   P     : B-spline patch
%   r     : max order of derivative
%   el    : left end if true, right end otherwise (optional)
%
% OUTPUT
%   K     : differentiation matrix at end point up r-th order

if nargin < 3
   el = true;
end
K = zeros(P.n, r+1);
if el
   K(1, 1) = 1;
   if P.p > 0 && r > 0
      r = min(P.p, r);
      for q = P.p-r+1:P.p
         j = 1:r-P.p+q;
         for i = r-P.p+q+1:-1:1
            l = P.U(i+P.p+1) - P.U(i+P.p+1-q);
            if l > 0
               K(i, j+1) = - q / l * K(i, j);
            end
            l = P.U(i+P.p) - P.U(i+P.p-q);
            if i > 0 && l > 0
               K(i, j+1) = K(i, j+1) + q / l * K(i-1, j);
            end
         end
      end
   end
else
   K(P.n, 1) = 1;
   if P.p > 0 && r > 0
      r = min(P.p, r);
      for q = P.p-r+1:P.p
         j = 1:r-P.p+q;
         for i = P.n-r+P.p-q:P.n
            l = P.U(i+q) - P.U(i);
            if l > 0
               K(i, j+1) = q / l * K(i, j);
            end
            l = P.U(i+q+1) - P.U(i+1);
            if i < P.n && l > 0
               K(i, j+1) = K(i, j+1) - q / l * K(i+1, j);
            end
         end
      end
   end
end

end
