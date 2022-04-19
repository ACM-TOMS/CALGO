function M = B_differentiation_all(P, r, xx, cl)

% Differentiation of all B-splines in given points
%
% INPUT
%   P     : B-spline patch
%   r     : order of derivative
%   xx    : vector of evaluation points
%   cl    : closed domain if true (optional)
%
% OUTPUT
%   M     : differentiation matrix 

if nargin < 4
   cl = true;
end
if r <= 0
   M = B_evaluation_all(P, xx, cl);
elseif r <= P.p
   tol = 1e-12;
   Pr = struct('p', P.p-r, 'n', P.n+r, 'U', P.U);
   Mr = B_evaluation_all(Pr, xx, cl);
   for q = P.p-r+1:P.p
      for i = P.p-q+1:P.n
         j = (xx >= P.U(i)-tol & xx <= P.U(i+q+1)+tol);
         l = P.U(i+q) - P.U(i);
         if l > 0 && any(j)
            Mr(i, j) = q / l * Mr(i, j);
         end
         l = P.U(i+q+1) - P.U(i+1);
         if l > 0 && any(j)
            Mr(i, j) = Mr(i, j) - q / l * Mr(i+1, j);
         end
      end
   end
   M = Mr(1:P.n, :);
else
   M = zeros(P.n, length(xx));
end

end
