function [ll,lld] = loglik(theta, X, p, q, r, miss)
  [A, B, Sig, mu] = theta2mat(theta, p, q, r);
  if nargin<6, miss = isnan(X); end
  if nargout==1
    [ll, ok] = varma_llm(X, A, B, Sig, mu, miss);
  else
    [ll, ok, lld] = varma_llm(X, A, B, Sig, mu, miss);
  end
  ascertain(ok);
end
