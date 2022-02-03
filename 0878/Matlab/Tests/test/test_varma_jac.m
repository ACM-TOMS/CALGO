%TEST_VARMA_JAC  Test Jacobian feature of varma_llc and varma_llm
%
%  TEST_VARMA_JAC compares gradients of varma_llc and varma_llm with numerical
%  gradient when change of variables is used. It also compares the likelihood
%  function value with the one obtained when no Jacobian is specified.

function test_varma_llc_jac
  fprintf('TESTING JACOBIAN FEATURE OF VARMA_LLC and VARMA_LLM...\n');
  ncase = testcase('number');
  randn('state',3); rand('state',3);
  for r=1:2:3
    for pq = {0 1; 0 2; 1 0; 1 2; 2 0; 2 1}'
      [p,q] = deal(pq{:});
      n = p+q+3;
      X = rand(r,n);
      m1 = r*r*p;      % A
      m2 = r*r*q;      % B
      m3 = r*(r+1)/2;  % Sig
      n1 = ceil(m1/5); % thA
      n2 = ceil(m2/5); % thB
      n3 = 2;          % thSig
      S1 = gallery('lehmer',r);
      S2 = gallery('minij',r);
      th = (rand(n1 + n2 + n3, 1))/(p+q)/r/r;
      Ja = rand(r^2*p, n1);
      Jb = rand(r^2*q, n2);
      Jc = [vech(S1) vech(S2)];
      J = blkdiag(Ja, Jb, Jc);
      Xm = makemissing(X, 8);
      thmu = [th; 0.1*(1:r)'];
      d(1) = diff_test(@fun_c, th, X, p, q, r, J);
      d(2) = diff_test(@fun_m, thmu, Xm, p, q, r, J);
      fprintf('  p,q,r = %d %d %d;  max-diff = %.1e\n',p,q,r,max(d));
      ascertain(max(d) < 1e-8);
    end
  end
  disp('  OK')
end

function [ll,lld] = fun_c(theta, X, p, q, r, J)
  [A, B, Sig] = theta2mat(theta, p, q, r, J);
  ll = varma_llc(X, A, B, Sig);
  if nargout>1
    [ll1, ok, lld] = varma_llc(X, A, B, Sig, J);
    ascertain(ok);
    ascertain(almostequal(ll, ll1));
  end
end

function [ll,lld] = fun_m(theta, X, p, q, r, J)
  miss = isnan(X);
  [A, B, Sig, mu] = theta2mat(theta, p, q, r, J);
  ll = varma_llm(X, A, B, Sig, mu, miss);
  if nargout>1
    [ll1, ok, lld] = varma_llm(X, A, B, Sig, mu, miss, J);
    ascertain(ok);
    ascertain(almostequal(ll, ll1));
  end
end

