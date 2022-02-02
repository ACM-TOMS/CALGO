%TEST_VARMA_LLC  Check that varma_llc calculates correct likelihood.
%
%  TEST_VARMA_LLC checks that varma_llc returns the same values as likelihood
%  calculated directly from S matrix. During development it was also ascertained
%  that varma_llc gives the same result as a previously published program for
%  complete data likelihood. Information about this comparison and the programs
%  used are contained in the subdirectory compare_as311 (see readme.txt file
%  there).

function test_varma_llc
  fprintf('TESTING VARMA_LLC...\n');
  ncase = testcase('number');
  randn('state',1);
  disp(['  Case           varma_llc        full      rel.diff  res.diff'])
  for j=1:ncase-2
    [A, B, Sig, name] = testcase(j);
    r = size(Sig,1);
    n = 10;
    X = varma_sim(A, B, Sig, n);
    [l, ok, res] = varma_llc(X, A, B, Sig, 'res');
    ascertain(ok)
    [lS, resS] = llresS(A, B, Sig, X);
    ldiff = abs(l-lS)/abs(l);
    rdiff = norm(res-resS);
    fmt = '  %-12s %13.9f %13.9f %9.1e %9.1e\n';
    fprintf(fmt, name, l, lS, rdiff, ldiff);
    ascertain(ldiff<1e-14 && rdiff<1e-13);
  end
  disp('  OK');
end

function [lS, resS] = llresS(A, B, Sig, X)
  % likelihood and residuals directly from S matrix
  n = size(X,2);
  r = size(Sig,1);
  q = length(B)/r;
  x = X(:);
  [C,G,W,S] = find_CGWS(A, B, Sig);
  SS = S_build(S, A, G, n);
  lS = (-n*r*log(2*pi) - log(det(SS)) - x'*(SS\x))/2;
  CC = cell(n,n);
  for k=1:n
    for j=1:n, CC{k,j} = zeros(r,r); end
    for j=0:min(q,k-1)
      CC{k-j, k} = C{j+1}';
    end
  end
  CC = cell2mat(CC);
  resS = reshape(CC*(SS\x), r, n);
end
