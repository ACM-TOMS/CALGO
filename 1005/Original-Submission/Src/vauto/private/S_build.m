%S_BUILD  Calculate covariance of x
%
%  SS = S_BUILD(S,A,G,n) finds the covariance matrix of the vector x =
%  [x1'...xn']' of all the values of a VARMA time series,
%
%                  x(t) = A1·x(t-1) + ... + Ap·x(t-p) + y(t)
%  where
%                  y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q),
%
%  and x(t), y(t) and eps(t) are r-dimensional with eps(t) N(0,Sig). S = 
%  {S0,S1...Sp}, with Sj = Cov(x(t),x(t-j)) and G = {G0 G1...Gq} with Gj =
%  Cov(y(t),x(t-j)) should be previously calculated with find_CGW. A should
%  contain [A1...Ap]. SS is the matrix:
%
%                   S0  S1' S2'...Sn-1'
%                   S1  S0  S1'...Sn-2'
%                   S2               :
%                   :                :
%                   Sn-1 ...... S1  S0
%
%  and the Sj are found with the recurrence relation:
%
%      Sj = A1*S(j-1) + A2*S(j-2) + ........ + Ap*S(j-p) + Gj
%
%  with Gj = 0 for j > q.

function SS = S_build(S, A, G, n)
  r = size(G{1},1);
  p = length(S) - 1;
  q = length(G) - 1;
  S = [S cell(1,n-p-1)];
  A = makecell(A);
  for j=p+1:n-1
    if j <= q, S{j+1} = G{j+1}; else S{j+1} = zeros(r,r); end
    for i=1:p
      S{j+1} = S{j+1} + A{i}*S{j+1-i};
    end
  end
  St = cell(1,n-1);
  Sc = cell(n-1,1);
  S0 = S{1};
  for i=1:n-1
    St{i} = S{i+1}';
    Sc{i} = S{i+1};
  end
  SS = cell(n,n);
  for i=1:n
    SS(i,i) = {S0};
    SS(i,i+1:end) = St(1:n-i);
    SS(i+1:end,i) = Sc(1:n-i);
  end
  SS = cell2mat(SS);
end
