%LAMBDA_T_MULTIPLY  Multiply vector with transpose of Lambda
%
%  Y = LAMBDA_T_MULTIPLY(A, X, r, n) returns Lam'·X in Y, where Lam is an
%  n·r×n·r lower block-band matrix given by:
%
%               Lam =   I
%                           I
%                               I
%                                   I
%                      -Ap ... -A2 -A1  I
%                          -Ap     -A2 -A1  I
%                              ...         ..  .. 
%                                  -Ap ...    -A1  I
%
%  and X is a vector with compatible dimensions. A should be {A1...Ap}. 
%
%  Y = LAMBDA_T_MULTIPLY(A, X, r, n, miss) multiplies with Lam(obs,obs)' where
%  obs are the observed indices in {1,2,...,r·n}. This call is used for
%  multiplying with the matrix Lam_o'.

function Y = lambda_T_multiply(A, X, r, n, miss)
  if nargin < 5, miss = false(r,n); end
  A = makecell(A);
  p = length(A);
  ko = find_missing_info(miss);
  Y = X;
  for t=p+1:n
    K = find(~miss(:,t));
    K1 = ko(t)+1 : ko(t+1);
    for k = 1:p
      L = find(~miss(:,t-k)); 
      L1 = ko(t-k)+1 : ko(t-k+1);
      Y(L1) = Y(L1) - A{k}(K,L)'*X(K1);
    end
  end
end
