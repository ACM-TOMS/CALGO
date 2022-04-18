% TEST_FIND_SM  Check that find_Sm works correctly.
%
%   TEST_FIND_SM asserts that find_Sm returns the right thing for testcases 1 to
%   7 (see testcase.m)
%
%   TEST_FIND_SM(k) checks only testcase k

function test_somsm_build(caseno)
  fprintf('TESTING FIND_SM...');
  if nargin<1, cases=1:7; else cases=caseno; end
  for tcase = cases
    [A, B, Sig, p, q, r, name] = testcase(tcase);
    n = p+q+3;
    [C, G, W] = find_CGW(A, B, Sig);
    LUvyw = vyw_factorize(A);
    S = vyw_solve(A, LUvyw, G);
    Scol = S_extend(A, G, S, n);
    SS = S_build(S, A, G, n);
    nPar = r^2*(p+q) + r*(r+1)/2;
    %
    % Call the function once (useful for debug):
    miss = false(r, n); miss([1 3]) = true;
    [f,g] = sommfun(mat2theta(A, B, Sig), p, q, r, n, miss);
    %
    % First try with nothing missing:
    Sm = find_Sm(Scol, false(r,n));
    ascertain(isempty(Sm));
    %
    % Check sizes of derivatives:
    [CCd, GGd, WWd] = find_CGW_deriv(A, B, Sig);
    RHS = vyw_deriv_rhs(A, GGd, S);
    Sd = vyw_solve(A, LUvyw, RHS);
    [G, Gd] = der2array(GGd);
    [Scol,Scold] = S_extend(A, G, S, n, Gd, Sd);
    miss = false(r,n); miss([1,3,4]) = true;
    [Sm, Smd] = find_Sm(Scol, miss, Scold);
    ascertain(isequal(size(Smd), [3, 3, nPar]));
    %
    % Now two cases with ca. 25% missing:
    for k = 1:2
      miss = false(r,n);
      miss(unique(floor(r*n*rand(round(r*n/4),1)) + 1)) = true;
      Sm = find_Sm(Scol, miss);
      ascertain(almostequal(Sm, SS(miss,miss)));
    end
    %
    % Finally compare derivatives with numerical ones:
    [d,g,gnum] = diff_test(@sommfun, mat2theta(A, B, Sig), p, q, r, n, miss);
    ascertain(d<1e-8);
  end
  disp('OK');
end

function [f,g] = sommfun(theta, p, q, r, n, miss);
  [A, B, Sig] = theta2mat(theta,p,q,r);
  nPar = r^2*(p+q) + r*(r+1)/2;
  if nargout == 1
    [C, G, W] = find_CGW(A, B, Sig);
    PLU = vyw_factorize(A);
    S = vyw_solve(A, PLU, G);
    Scol = S_extend(A, G, S, n);
    Sm = find_Sm(Scol, miss);
    f = Sm(:);
  else
    nmiss = sum(miss(:));
    nobs = n*r - nmiss;
    [C, G, W] = find_CGW(A, B, Sig);
    [CCd, GGd, WWd] = find_CGW_deriv(A, B, Sig);
    PLU = vyw_factorize(A);
    S = vyw_solve(A, PLU, G);
    RHS = vyw_deriv_rhs(A, GGd, S);
    Sd = vyw_solve(A, PLU, RHS);
    [G, Gd] = der2array(GGd);
    [Scol,Scold] = S_extend(A, G, S, n, Gd, Sd);
    [Sm, Smd] = find_Sm(Scol, miss, Scold);
    f = Sm(:);
    g = reshape(Smd, nmiss^2, []);
  end
end  
