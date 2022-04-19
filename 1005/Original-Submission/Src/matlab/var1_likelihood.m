%  VAR1_LIKELIHOOD  Calculate the likelihood for a vector-autoregressive time
%  series model with one term (VAR(p) with p=1).
%
%    l = VAR1_LIKELIHOOD(X, A, Sig) returns the negative log-likelihood
%    [l, Aa, La] = VAR1_LIKELIHOOD(X, A, Sig) returns also the derivatives of l
%    with respect to A and the lower triangle of Sig.

function [l, Aa, Siga] = var1_likelihood(X, A, Sig)
  % COMPUTE FUNCTION VALUE
  [r,n] = size(X);
  Sig = symm(Sig);
  U = dlyap1(A, Sig);
  S = tril(U);
  LSig = chol(Sig, 'lower');          % potrf
  LS = chol(S, 'lower');              % potrf
  W = zeros(size(X));
  W(:,1) = X(:,1);
  W(:,2:n) = X(:,2:n) - A*X(:,1:n-1); % gemm
  W(:,1) = LS\W(:,1);                 % trsv
  W(:,2:n) = LSig\W(:,2:n);           % trsm
  w = W(:);
  v = dot(w,w);
  sumdiag = 0;
  for i=1:r
    sumdiag = sumdiag + log(LS(i,i)) + (n-1)*log(LSig(i,i)); 
  end
  l = -(n*r*log(2*pi) + 2*sumdiag + v)/2;
  if nargout == 1, return, end
  
  % COMPUTE ADJOINT
  la = 1;
  va = -la/2;
  sumdiaga = -la;
  LSa = zeros(r,r);
  LSiga = zeros(r,r);
  for i=r:-1:1
    LSiga(i,i) = (n-1)*sumdiaga/LSig(i,i);
    LSa(i,i) = sumdiaga/LS(i,i);
  end
  wa = 2*va*w;                                 % dot_rmd
  Wa = reshape(wa, r, n);
  Wa(:,2:n) = LSig'\Wa(:,2:n);                 % trsm_rmd
  LSiga = LSiga - tril(Wa(:,2:n)*W(:,2:end)'); % trsm_rmd
  Wa(:,1)   = LS'\Wa(:,1);                     % trsv_rmd
  LSa = LSa - tril(Wa(:,1)*W(:,1)');           % trsv_rmd
  Aa = -Wa(:,2:n)*X(:,1:n-1)';                 % gemm_rmd
  Sa = chol_rmd(S, LS, zeros(r), LSa);         % potrf_rmd
  Siga = chol_rmd(Sig, LSig, zeros(r), LSiga); % potrf_rmd
  Ua = (symm(Sa) + dg(Sa))/2;                   % isym_rmd
  V = dlyap1(A', Ua);                          
  Siga = Siga + 2*tril(V) - dg(V);             % sym_rmd
  Aa = Aa + 2*V*A*U;
end
