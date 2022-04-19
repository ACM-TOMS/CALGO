%TEST_VARMA_LLC Compare varma_llc with as311
%
%  TEST_VARMA_LLC checks that varma_llc returns the same value as as311 for all
%  test cases that "testcase" creates. Also compares with likelihood calculated
%  directly from S matrix.

function test_varma_llc
  disp('COMPARING VARMA_LLC WITH ALGORITHM AS-311');
  ncase = testcase('number');
  randn('state',1);
  disp('  Case            varma_llc       as311     rel.diff')
  for j=1:ncase-2
    [A, B, Sig, name] = testcase(j);
    r = size(Sig,1);
    n = 10;
    X = varma_sim(A, B, Sig, n);
    [logelf,F1,F2,a,ifault,Sig] = as311(X, A, -B, Sig, [], 1, 1, -1);
    l = varma_llc(X, A, B, Sig, 'res');
    ldiff = abs(l-logelf)/abs(l);
    fmt = '  %-12s %13.9f %13.9f %9.1e\n';
    fprintf(fmt, name, l, logelf, ldiff);
  end
end
