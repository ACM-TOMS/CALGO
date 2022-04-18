% TEST_FIND_LAMBDA_OM  Check that lambda_om is calculated correctly
%
%   TEST_FIND_LAMBDA_OM checks that functions find_Acol and find_lambda_om
%   return the correct values. Runs testcases 1-9
%
%   TEST_FIND_LAMBDA_OM(k) tries only testcase k

function test_find_laom(caseno)
  fprintf('TESTING FIND_ACOL AND FIND_LAMBDA_OM... ');
  randn('state',1); rand('state',1);
  if nargin<1, cases=1:9; else cases=caseno; end
  for tcase = cases
    [A, B, Sig, p, q, r, name] = testcase(tcase);
    n = p+q+3;
    Acol = find_Acol(A, r);
    Lam = lambda_build(A, r, n);
    %
    % Nothing missing:
    miss = false(r,n);
    Laom = find_lambda_om(Acol, miss);
    ascertain(isempty(Laom));
    %
    % Two cases with missing cases (fraction missing given by PCT):
    PCT = [0.10 0.25];
    for k = 1:2
      miss(unique(floor(r*n*rand(round(r*n*PCT(k)),1)) + 1)) = true;
      obs = ~miss;
      Laom = find_lambda_om(Acol, miss);
      LaomOK = Lam(obs,miss);
      mLaom = size(Laom, 1);
      ascertain(equal(Laom, LaomOK(1:mLaom, :)));
      ascertain(all(all(LaomOK(mLaom+1:end,:)==0)))
    end
    %
    % Compare derivatives with numerical ones:
    [d,g,gnum] = diff_test(@laomfun, mat2theta(A, B, Sig), p, q, r, miss);
    % g=g{1},gnum=gnum{1}
    ascertain(d<1e-8);
  end
  disp('OK');
end

function [f,g] = laomfun(theta, p, q, r, miss);
  [A, B, Sig] = theta2mat(theta,p,q,r);
  nPar = r^2*(p+q) + r*(r+1)/2;
  if nargout == 1
    Acol = find_Acol(A, r);
    Laom = find_lambda_om(Acol, miss);
  else
    [Acol, Acold] = find_Acol(A, r, nPar);
    [Laom, Laomd] = find_lambda_om(Acol, miss, Acold);
    g = reshape(Laomd,[],nPar);
  end
  f = Laom(:);
end

function e = equal(x, y)
  e = isequal(size(x), size(y)) && isequal(x(:), y(:));
end
