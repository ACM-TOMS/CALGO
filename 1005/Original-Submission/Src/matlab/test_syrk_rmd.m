% TEST_SYRK_RMD  Check that syrk_rmd works as documented.

function test_syrk_rmd
  for m=1:5
    for n=1:5
      test(m, n);
    end
  end
end

function test(m, n)
  % Test adjoint of mxn matrix times n vector with nf outputs
  A = randPM(m, n); 
  u = vec(A);
  ftest(@oper, @rmd, u, 'repeats', 2, 'testcase', {m, n, 'N'});
  ftest(@oper, @rmd, u, 'repeats', 2, 'testcase', {m, n, 'T'});
end

function v = oper(u, m, n, trans) % syrk operation
  trsp = nargin > 3 && strcmpi(trans, 't');
  A = reshape(u, m, n);
  if trsp
    L = tril(A'*A);
  else
    L = tril(A*A');
  end
  v = vech(L);
end

function ua = rmd(u, ~, ua, va, m, n, trans)
  A = reshape(u, m, n);
  Aa = reshape(ua, m, n);
  La = ivech(va);
  Aa = syrk_rmd(A, Aa, La, trans);
  ua = vec(Aa);
end
