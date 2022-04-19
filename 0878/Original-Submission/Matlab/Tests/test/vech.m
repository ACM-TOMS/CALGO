%VECH  Change lower triangle to column vector
%
%  v = VECH(A), where A is an n×n matrix returns an n·(n+1)/2 dimensional
%  column vector with the columns of the lower triangle of A placed one after
%  another.
%
%  V = VECH(A) where A is n×n×N returns an n·(n+1)/2 × N dimensional matrix
%  with VECH(A(:,:,j)) in its j-th column.

function v = vech(A)
  if isempty(A), v=[];
  else
    [n,m,N] = size(A);
    ascertain(n==m);
    v = zeros(n*(n+1)/2, N);
    m = 1;
    for i=1:n
      m1 = m + n-i;
      v(m:m1, :) = A(i:n, i, :);
      m = m1 + 1;
    end
  end
end
