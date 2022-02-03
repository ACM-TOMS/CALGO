function ss = MDB_evaluation_spline(MP, H, cc, xx)

% Evaluation of a multi-degree spline in given points
%
% INPUT
%   MP    : MDB-spline multi-patch
%   H     : extraction matrix
%   cc    : vector of coefficients
%   xx    : vector of evaluation points
%
% OUTPUT
%   ss    : vector of spline evaluation values

yy = xx;
m = length(MP.P) - 1;
dd = reshape(cc, 1, []) * H;
ss = zeros(1, length(xx));
for i = 1:m
   dd_loc = dd(MP.mu(i)+1:MP.mu(i+1));
   j = (yy >= MP.P(i).U(1) & yy <= MP.P(i).U(end));
   ss(j) = B_evaluation_spline(MP.P(i), dd_loc, yy(j));
   yy = yy - MP.P(i).U(end) + MP.P(i+1).U(1);
end
dd_loc = dd(MP.mu(m+1)+1:MP.mu(end));
j = (yy >= MP.P(end).U(1) & yy <= MP.P(end).U(end));
ss(j) = B_evaluation_spline(MP.P(end), dd_loc, yy(j));

end
