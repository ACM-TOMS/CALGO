function M = MDB_evaluation_all(MP, H, xx)

% Evaluation of all MDB-splines in given points
%
% INPUT
%   MP    : MDB-spline multi-patch
%   H     : extraction matrix
%   xx    : vector of evaluation points
%
% OUTPUT
%   M     : evaluation matrix 

yy = xx;
m = length(MP.P) - 1;
N = zeros(MP.mu(end), length(xx));
for i = 1:m
   N(MP.mu(i)+1:MP.mu(i+1), :) = B_evaluation_all(MP.P(i), yy, false);
   yy = yy - MP.P(i).U(end) + MP.P(i+1).U(1);
end
N(MP.mu(m+1)+1:MP.mu(end), :) = B_evaluation_all(MP.P(end), yy, true);
M = H * N;

end
