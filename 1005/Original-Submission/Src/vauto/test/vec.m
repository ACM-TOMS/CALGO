%VEC  Change matrix to column vector
%
%  v = VEC(A), where A is an m × n matrix, returns the columns of A one after
%  another. It is equivalent to v=A(:).
%
%  v = VEC(A), where A is m × n × N, returns an m·n × N matrix with
%  vec(A(:,:,j)) in its j-th column.

function v = vec(A)
  [m,n,N] = size(A);
  if N==1
    v = A(:);
  else
    v = zeros(m*n, N);
    for j=1:N, v(:,j) = reshape(A(:,:,j), m*n, 1); end
  end
end 
