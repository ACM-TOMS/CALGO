%TEST_VARMA_LLM_DERIV
%
%  TEST_VARMA_LLM_DERIV compares gradient of missing value log-likelihood
%  function with a numerical gradient for a few testcases.

function test_varma_llm_deriv
  disp('TESTING VARMA_LLM_DERIV...')
  disp('  Max analytic / numerical gradient difference:')
  rand('state',2); randn('state',3);
  for j = 1:8
    [A, B, Sig, p, q, r, name] = testcase(j);
    mu = 0.1*(1:r)';
    n = p+q+4;
    X = varma_sim(A, B, Sig, n, mu);
    if r>=3, I=[-1 0]; elseif r==2, I=[-4 0]; else I=[-4:0 20]; end
    for i=I
      Xm = makemissing(X,i);
      miss = isnan(Xm);
      theta = mat2theta(A, B, Sig, mu);
      d(i+5) = diff_test(@loglik, theta, Xm, p, q, r, miss);
      if d(i+5)>=1e-8, d(i+5), end
      ascertain(d(i+5)<1e-8)
    end
    fprintf('  Testcase %-12s %.1e\n', name, max(d));
  end
  disp('  OK')
end
