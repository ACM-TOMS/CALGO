%C_MULTIPLY  Multiply vector by expanded C-matrix
%
%  Y = C_MULTIPLY(C, X, n, miss) returns C*X in Y, where C is obtained from the
%  n·r×n·r lower block-band matrix:
%
%                       C0  C1' C2' ... Cq'
%                           C0  C1' C2' ... Cq'
%                               ...           ...
%                                   C0  C1' ... Cq'
%                                       ...      :
%                                           C0  C1'
%                                               C0
%
%  by removing columns with indices in miss (which is a subset of {1,2,...,r·n}.
%  C = {C0 C1...Cq} is as returned by find_CGW (Cj is r×r).

function Y = C_multiply(C, X, n, miss)
  q = length(C) - 1;
  r = size(C{1}, 1);
  if nargin < 4, miss = false(r,n); end
  ko = find_missing_info(miss);
  Y = zeros(r, n);
  for t = 1:n
    K = find(~miss(:,t));
    K1 = (ko(t)+1 : ko(t+1))';
    for k = min(q,t-1):-1:0
      Y(:,t-k) = Y(:,t-k) + C{k+1}(K,:)'*X(K1);
    end
  end
  Y = reshape(Y, n*r, 1);
end
