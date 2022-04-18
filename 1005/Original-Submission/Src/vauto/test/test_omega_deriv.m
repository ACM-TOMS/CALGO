% TEST_OMEGA_DERIV  Test complete data derivatives of Omega and its factors
%
%   TEST_OMEGA_DERIV checks the functions omega_build_deriv, omega_factor_deriv,
%   omega_forward_deriv in the complete data case. In addition the derivative
%   features of lambda_multiply, omega_ltl, omega_back_sub and omega_logdet are
%   tested, as well as the derivative of z'z (also only for complete data). Thus
%   all the ingredients in the derivative of the log-likelihood function are
%   tested, thus providing a double-check in addition to the checks of
%   test_varma_llc_deriv and test_varma_llm_deriv.
%
%   By default a set of mostly r=2 testcases are generated and used. Specify
%   TEST_OMEGA_DERIV(n) to use the n-th named case from "testcase".
%
%   Use "TEST_OMEGA_DERIV quiet" to suppress output.

function test_omega_deriv(varargin)
  [varargin,quiet] = getflags(varargin,'quiet');
  randn('state', 8); rand('state', 8);
  pj = [0 1 1 1 1 2 2 2 2 2];
  qj = [2 0 0 1 3 0 1 2 3 0];
  rj = [1 1 2 2 2 2 2 2 2 3];
  fprintf('TESTING DERIVATIVES OF OMEGA, LL'', L''L, w, z, AND z''z...');
  fprintf_if(~quiet,'\n');
  fmt = '    %-28s = %.1e\n';
  ncase = testcase('number');
  if ~isempty(varargin) 
    cases = varargin; 
  else
    cases = num2cell([pj;qj;rj]);
  end
  for tcase = cases
    [A, B, Sig, name] = testcase(tcase{:});
    theta = mat2theta(A, B, Sig);
    [p,q,r] = get_dimensions(A, B, Sig);
    n = p+q+2;
    fprintf_if(~quiet, '  Testcase p=%d q=%d r=%d:\n', p, q, r);
    [dmax,g,gnum] = diff_test(@fun, theta, p, q, r, n);
    fprintf_if(~quiet, fmt, 'Chol-factor max rel.diff', dmax(1));
    fprintf_if(~quiet, fmt, 'ltl max rel.diff', dmax(2));
    fprintf_if(~quiet, fmt, 'lambda_multiply max rel.diff', dmax(3));
    fprintf_if(~quiet, fmt, 'forward max rel.diff', dmax(4));
    fprintf_if(~quiet, fmt, 'backsub max rel.diff', dmax(5));
    fprintf_if(~quiet, fmt, 'logdet max rel.diff', dmax(6));
    fprintf_if(~quiet, fmt, 'z''z max rel.diff', dmax(7));
    %gnum{2},g{2}
    ascertain(all(dmax < 2e-8));
  end
  disp('  OK');
end

function [f,g] = fun(theta,p,q,r,n)
  nPar = r^2*(p+q) + r*(r+1)/2;
  ko = 0:r:r*n;
  [A, B, Sig] = theta2mat(theta, p, q, r);
  [C, G, W, S] = find_CGWS(A, B, Sig);
  [Su, Olow] = omega_build(S, G, W, p, n);
  [Lu, Ll, info] = omega_factor(Su, Olow, p, q, ko);
  if info ~= 0, error('Omega not positive definite'); end
  [Mu, Ml, info] = omega_ltl(Su, Olow, p, q, ko); ascertain(info==0);
  X = reshape(1:r*n, r, n);
  w = lambda_multiply(A, X(:), false(r, n));
  z = omega_forward(Lu, Ll, w, p, q, ko);
  zb = omega_back_sub(Lu, Ll, w, p, q, ko);
  ld = omega_logdet(Lu, Ll, p, q, ko);
  f{1} = [vech(Lu); getvec(Ll, r)];
  f{2} = [vech(Mu); getvec(Ml, r)];
  f{3} = w;
  f{4} = z;
  f{5} = zb;
  f{6} = ld;
  f{7} = z'*z;
  if nargout>1
    [C, G, W, S, Cd, Gd, Wd, Sd] = find_CGWS(A, B, Sig);
    [Sud, Olowd] = omega_build_deriv(Sd, Gd, Wd, p, n);
    [Lud, Lld] = omega_factor_deriv(Sud, Olowd, Lu, Ll, p, q, ko);
    [Mu, Ml, info, Mud, Mld] = omega_ltl(Su, Olow, p, q, ko, Sud, Olowd);
    xd = zeros(r*n, 1, nPar);
    [w, wd] = lambda_multiply(A, X(:), false(r, n), xd);
    zd = omega_forward_deriv(Lu, Ll, Lud, Lld, z, wd, p, q, ko);
    [zb,zbd] = omega_back_sub(Lu, Ll, w, p, q, ko, Lud, Lld, wd);
    [ld,ldd] = omega_logdet(Lu, Ll, p, q, 0:r:r*n, Lud, Lld);
    g{1} = [vech(Lud); getvec(Lld,r)];
    g{2} = [vech(Mud); getvec(Mld,r)];
    g{3} = [squeeze(wd)];
    g{4} = squeeze(zd);
    g{5} = squeeze(zbd);
    g{6} = ldd;
    g{7} = 2*z'*squeeze(zd);
  end
end

function x = getvec(Olow, r)
  % Convert Olow matrix to a vector.
  m = size(Olow,1)/r;
  q = size(Olow,2)/r - 1;
  M = size(Olow,3);
  N1 = r^2*q;
  N2 = N1 + r*(r+1)/2;
  x = zeros(N2, m, M);
  K = 1:r;
  for k = 1:m
    if q>0, x(1:N1, k, :) = reshape(Olow(K, 1:q*r, :), N1,1,[]); end
    x(N1+1:end, k, :) = vech(Olow(K, q*r+1:end, :));
  end
  x = reshape(x, N2*m, M);
end
