%TEST_VYW_DERIV Test derivative of vyw solutions
%  TEST_VYW_DERIV checks the derivative returned by vyw_deriv_rhs and subsequent
%  vyw_solve (for the derivative of S) against numerical differentiation.
%  TEST_VYW_DERIV QUIET minimizes output and TEST_VYW_DERIV(k) runs testcase k.
%

function test_vyw_deriv(varargin)
  [varargin, quiet] = getflags(varargin,'quiet');
  fprintf('TESTING VYW_DERIV...'); fprintf_if(~quiet,'\n');
  tcspec = ~isempty(varargin);
  if tcspec
    cases = varargin{1};
  else
    pj = [1 1 2 2 2 2 3 2];
    qj = [0 1 0 1 2 3 4 0];
    rj = [1 2 2 2 2 2 2 3];
    cases=1:length(pj);
  end
  for j=cases
    if tcspec, [A, B, Sig, p, q, r, name] = testcase(j);
    else       [A, B, Sig, p, q, r, name] = testcase(pj(j),qj(j),rj(j)); end
    theta = mat2theta(A, B, Sig);
    fprintf_if(~quiet, '  Testcase p=%d q=%d r=%d: ', p, q, r);
    dmax = diff_test(@fun, theta, p, q, r);
    fprintf_if(~quiet, 'max rel.diff=%.1e\n', dmax);
    ascertain(dmax<1e-8);
  end
  disp('  OK');

  function [f,g]=fun(theta,p,q,r)
    N = length(theta);
    [A, B, Sig] = theta2mat(theta,p,q,r);
    [C, G, W] = find_CGW(A, B, Sig);
    PLU = vyw_factorize(A);    
    S = vyw_solve(A, PLU, G);
    f0 = vech(S{1});
    f1 = reshape(cell2mat(S(2:end-1)), r^2*(p-1), 1);
    f = [f0; f1];
    if nargout>=2
      [CCd, GGd, WWd] = find_CGW_deriv(A, B, Sig);
      RHS = vyw_deriv_rhs(A, GGd, S);
      Sd = vyw_solve(A, PLU, RHS);
      Sd0 = vech(Sd{1});
      Sd1 = reshape(cell2mat(Sd(2:end-1)), r^2*(p-1), N);
      g = [Sd0; Sd1];
    end
  end
end
