function test_var1_likelihood(r, n)
  if nargin == 0
    for r = 1:4
      for n = 1:4
        test_var1_likelihood(r,n)
      end
    end
  end
  A = randPM(r, r)/r^2;
  Sig = gallery('lehmer', r);
  X = randPM(r, n);
  u = [vec(A); vech(Sig)];
  accum = false(size(u));
  ftest(@oper, @rmd, u, 'accum', accum, 'tol', 1e-6, 'testcase', {r, X});
end

function v = oper(u, n, X)
  A = reshape(u(1:n^2), n, n);
  L = ivech(u(n^2+1:end));
  v = var1_likelihood(X, A, L);
end

function ua = rmd(u, ~, ~, ~, n, X)
  A = reshape(u(1:n^2), n, n);
  L = ivech(u(n^2+1:end));
  [~, Aa, La] = var1_likelihood(X, A, L);
  ua = [vec(Aa); vech(La)];
end
