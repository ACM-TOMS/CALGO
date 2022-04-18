% DLYAP1  Solve discrete Lyapunov equation
% 
% S = DLYAP1(A, Sig) solves the discrete Lyapunov equation S - A*S*A' = Sig. for
% symmetric Sig. It uses Gaussian elimination to solve the equation directly,
% and is adapted from  vyw_factorize/vyw_solve in [1]. It has the advantages
% that it is not GPL-licensed and that it is easy to translate to Fortran.
% However it is considerably slower than the Matlab control system toolbox dlyap
% (4 times for r = 10, 100 times for r = 100).
%
% [1] K Jonasson: Algorithm 878: Exact VARMA likelihood and its gradient
%     for complete and incomplete data with Matlab, ACM TOMS 2008.
%
% Copyright 2017. Kristján Jónasson. Free to use as per the MIT license.  

function S = dlyap1(A, Sig)
  n = size(A,1);
  nn = n*(n+1)/2;
  F = zeros(nn,nn);
  sig = zeros(nn,1);
  j1 = 1;
  j2 = n;
  for j = 1:n
    i1 = 1;
    i2 = n;
    Aj = A(:,j:end);
    for i = 1:n
      F(i1:i2, j1:j2) = - A(i,j)*Aj(i:end,:) - A(i:end,j)*Aj(i,:);
      i1 = i2 + 1;
      i2 = i2 + n - i;
    end
    F(j1,j1) = F(j1,j1) + 1;
    sig(j1:j2) = Sig(j:end,j);
    j1 = j2 + 1;
    j2 = j2 + n - j;
  end
  F = F + eye(nn);
  %[L,U,p] = lu(F,'vector');
  %PLU = {p,L,U};
  s = F\sig;
  S = zeros(n,n);
  j1 = 1;
  j2 = n;
  for i=1:n
    S(i:end,i) = s(j1:j2);
    j1 = j2 + 1;
    j2 = j2 + n - i;
  end
  S = S + S';
end
