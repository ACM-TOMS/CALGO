%TEST_VAR_LL  Test the pure autoregressive likelihood function var_ll
%
%  TEST_VAR_LL compares likelihood calculated with var_ll with that calculated
%  with varma_llm for test cases that "testcase" creates. The gradient is also
%  compared.

function test_var_ll(varargin)
  fprintf('TESTING VAR_LL...');
  if ~isempty(varargin)
    cases = {varargin(1)};
    if length(varargin) > 1, patt = varargin{2}; end
  else
    cases = num2cell(num2cell(1:testcase('number')-1));
  end
  ncase = testcase('number');
  randn('state',3); rand('state',3);
  maxdiff = 0;
  % Check non-stationary models:
  A = [0.5 0.51]; B = []; 
  X = zeros(1,10); miss = isnan(X);
  mu = 0;
  Sig = -1; [l, ok] = var_ll(X, A, Sig, mu, miss); ascertain(~ok);  
  Sig = 1;  [l, ok] = var_ll(X, A, Sig, mu, miss); ascertain(~ok);
  for j=1:ncase-1
    [A, B, Sig, p, q, r, name] = testcase(j);
    if q==0
      mu = 0.01*(1:r)';
      r = size(Sig,1);
      n = p+3;
      X = varma_sim(A, [], Sig, n, mu);
      compare(A, Sig, X, p, r);
      mubar = repmat(mu, n, 1);
      for i=[20]
        X1 = makemissing(X,i);
        miss = isnan(X1);
        [l1, ok] = varma_llm(X1, A, [], Sig, mu, miss);
        [l, ok] = var_ll(X1, A, Sig, mu, miss); ascertain(ok)
        ascertain(almostequal(l, l1));
        [l2, ok, res, xm] = var_ll(X1, A, Sig, mu, miss, 'res_miss');
        ascertain(isequal(l, l2));
        [l2, ok, res1, xm1] = varma_llm(X1, A, [], Sig, mu, miss, 'res_miss'); 
        ascertain(ok);
        ascertain(almostequal(res1, res))
        ascertain(almostequal(xm1, xm));
        [l, ok, ld] = var_ll(X1, A, Sig, mu, miss);
        [l, ok, ld] = varma_llm(X1, A, [], Sig, mu, miss);
      end;
    end
  end
  disp('  OK');
end

function compare(A, Sig, X, p, r)
  % Compare function values of (a) var_ll with zero mu and all-false miss, (b)
  % var_ll without mu- and miss-parameters and (c) varma_llc. Also check
  % redidual calculation.
  miss = false(size(X));
  ll1 = var_ll(X, A, Sig, zeros(r,1), miss);
  ll2 = var_ll(X, A, Sig);
  ll3 = varma_llc(X, A, [], Sig);
  ascertain(ll1==ll2 && almostequal(ll3, ll1));
  [ll4, ok, res1] = var_ll(X, A, Sig, zeros(r,1), miss, 'res');
  [ll5, ok, res2] = var_ll(X, A, Sig, 'res');
  [ll6, ok, res3] = varma_llm(X, A, [], Sig, zeros(r,1), miss, 'res');
  ascertain(almostequal(ll4,ll1) && almostequal(ll5,ll2) && almostequal(ll6,ll3));
  ascertain(almostequal(res1, res2));
  ascertain(almostequal(res2, res3));
end
