%TEST_VARMA_LLM  Test the missing value likelihood function varma_llm
%
%  TEST_VARMA_LLM compares likelihood calculated with varma_llm with a direct
%                 calculation for all test cases that "testcase" creates and
%                 six different missing value patterns (from "makemissing").
%  TEST_VARMA_LLM QUIET runs quietly
%  TEST_VARMA_LLM(N) runs only the N-th named testcase from "testcase"
%  TEST_VARMA_LLM(N,M) runs N-th case on M-th missing value pattern
%  TEST_VARMA_LLM FULL runs a thorough test with more combinations of (p,r,q)
%                      and more missing value patterns
%

function test_varma_llm(varargin)
  [args,quiet] = getflags(varargin,'quiet');
  [args,full] = getflags(args,'full');
  fprintf('TESTING VARMA_LLM'); if full, fprintf(' (COMPREHENSIVE CHECK)'); end
  disp ' '
  patt = [];
  if full
    cases = {};
    for p=0:3, for q=0:3, for r=1:3, cases{end+1} = {p,q,r}; end, end, end
  elseif ~isempty(varargin)
    cases = {varargin(1)};
    if length(varargin) > 1, patt = varargin{2}; end
  else
    cases = num2cell(num2cell(1:testcase('number')-1));
  end
  randn('state',1);
  rand('state',1);
  maxdiff = 0;
  % Check non-stationary models:
  A = [0.5 0.51]; B = [];
  X = zeros(1,10); miss = isnan(X);
  mu = 0;
  Sig = -1; [l, ok] = varma_llm(X, A(1), B, Sig, mu, miss); ascertain(~ok);
  Sig = 1;  [l, ok] = varma_llm(X, A, B, Sig, mu, miss); ascertain(~ok);
  fmt = '  Testcase %s, max diff.=%.1e, max res.diff.=%.1e\n';
  for j = 1:length(cases)
    tcase = cases{j};
    if full, tcstring = sprintf('p=%d q=%d r=%d', tcase{:});
    else tcstring = int2str(tcase{1}); end
    [A, B, Sig, p, q, r, name] = testcase(tcase{:});
    mu = 0.01*(1:r)';
    r = size(Sig,1);
    n = max(p,q)+6;
    X = varma_sim(A, B, Sig, n, mu);
    mubar = repmat(mu, n, 1);
    [SS, CC] = SC_build(A, B, Sig, n);
    maxdiff(j) = 0; rdiff(j) = 0; mdiff(j) = 0;
    if ~isempty(patt), misspat = patt;
    elseif ~full,      misspat = [0 -4 -3 -2 -1 20];
    elseif  full,      misspat = [0 -4 -3 -2 -1 2 4 10 15 20 25]; end
    for i = misspat
      X1 = makemissing(X,i);
      miss = isnan(X1);
      if all(miss(:)==0) % CHECK THAT LLM AND LLC RETURN THE SAME:
        zmu = zeros(r,1);
        [l, ok] = varma_llm(X1, A, B, Sig, zmu, miss);
        ascertain(ok)
        l1 = varma_llc(X1, A, B, Sig);
        ascertain(almostequal(l, l1));
      end
      [l, ok] = varma_llm(X1, A, B, Sig, mu, miss);  ascertain(ok)
      ascertain(~isnan(l));
      % CHECK FOR SAME LIKELIHOOD WHEN RES AND/OR XM ARE ALSO RETURNED:
      [l1, ok, res, xm] = varma_llm(X1, A, B, Sig, mu, miss, 'res_miss');
      ascertain(ok);
      ascertain(almostequal(l, l1));
      if j==2 % check that only 'res' works
        [l, ok, res1] = varma_llm(X1, A, B, Sig, mu, miss, 'res');
        ascertain(ok && isequal(res,res1));
      end
      obs = ~miss;
      nObs = sum(obs(:));
      SSo = SS(obs,obs);
      xo = reshape(X(obs),[],1) - mubar(obs);
      lS = (-nObs*log(2*pi) - log(det(SSo)) - xo'*(SSo\xo))/2;
      resS = reshape(CC(:, obs)*(SSo\xo), r, n);
      xmS = SS(miss,obs)*(SSo\xo) + mubar(miss);
      maxdiff(j) = max(maxdiff(j), reldiff(l,lS));
      rdiff(j) = max(rdiff(j), reldiff(res,resS));
      if ~isempty(miss), mdiff(j) = max(mdiff(j), reldiff(xmS,xm)); end
    end;
    fprintf_if(~quiet, fmt, tcstring, maxdiff(j), rdiff(j));
  end
  fmt = '  Overall,  max diff.=%.1e, max res.diff.=%.1e, max xmiss diff=%.1e, ';
  fprintf(fmt, max(maxdiff), max(rdiff), max(mdiff));
  disp ' '
  disp '  ("max diff" is the difference between varma_llm function value and a'
  disp '  value calculated directly from the large covariance matrix S=Cov(x))'
  ascertain(max(maxdiff)<1e-14);
  ascertain(max(rdiff)<2e-14);
  ascertain(max(mdiff)<5e-14);
  disp('  OK');
end

function r = reldiff(a,b)
  ascertain(isequal(size(a),size(b)));
  if isempty(a)
    r = 0;
  else
    rmx = max([1,max(abs(a(:))), max(abs(b(:)))]);
    r = max(abs(a(:)-b(:)))/rmx;
  end
end

function [SS, CC] = SC_build(A, B, Sig, n) % cov(x,x) and cov(x,eps)
  [C,G,W,S] = find_CGWS(A, B, Sig);
  r = size(Sig,1);
  q = length(B)/r;
  SS = S_build(S, A, G, n);
  CC = cell(n,n);
  for k=1:n
    for j=1:n, CC{k,j} = zeros(r,r); end
    for j=0:min(q,k-1)
      CC{k-j, k} = C{j+1}';
    end
  end
  CC = cell2mat(CC);
end
