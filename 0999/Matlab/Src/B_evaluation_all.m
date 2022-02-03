function M = B_evaluation_all(P, xx, cl)

% Evaluation of all B-splines in given points
%
% INPUT
%   P     : B-spline patch
%   xx    : vector of evaluation points
%   cl    : closed domain if true (optional)
%
% OUTPUT
%   M     : evaluation matrix 

if nargin < 3
   cl = true;
end
tol = 1e-12;
M = zeros(P.n, length(xx));
for i = P.p+1:P.n
   j = (xx >= P.U(i) & xx < P.U(i+1));
   M(i, j) = 1;
   if cl && P.U(i+1) == P.U(end)
      M(i, abs(xx-P.U(end))<tol) = 1;
   end
end
for q = 1:P.p
   for i = P.p-q+1:P.n
      j = (xx >= P.U(i)-tol & xx <= P.U(i+q+1)+tol);
      l = P.U(i+q) - P.U(i);
      if l > 0 && any(j)
         M(i, j) = (xx(j) - P.U(i)) / l .* M(i, j);
      end
      l = P.U(i+q+1) - P.U(i+1);
      if l > 0 && any(j)
         M(i, j) = M(i, j) + (P.U(i+q+1) - xx(j)) / l .* M(i+1, j);
      end
   end
end

end
