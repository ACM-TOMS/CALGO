%TEST_LAMBDA_MULTIPLY  Tests function lambda_multiply
%
%  TEST_LAMBDA_MULTIPLY checks that lambda_multiply returns the right value for
%  a few test-cases. Further checks of lambda_multiply are made by
%  test_omega_deriv and test_var_ll_jac.
%
function test_lambda_multiply(varargin)
  [args, qui] = getflags(varargin,'quiet');
  fprintf('TESTING LAMBDA_MULTIPLY... '); 
  fprintf_if(~qui,'\n');
  r = 1;
  n = 4;
  x = (1:n)';
  A = 1;
  Lam = [1 0 0 0;-A 1 0 0;0 -A 1 0;0 0 -A 1];
  W = lambda_multiply([A], x, false(r, n));
  ascertain(almostequal(W, Lam*x));
  %
  % TRY A MISSING VALUE:
  obs = [1,2,4]; miss = true(r, n); miss(obs) = false;
  Wo = lambda_multiply(A, x(obs), miss);
  ascertain(almostequal(Wo, Lam(obs,obs)*x(obs)));
  x = (1:12)';
  n = 6; r = 2;
  O = zeros(2);
  A1 = [0 1;1 1];
  A2 = [2 3;2 2];
  A3 = [4 5;3 3];
  Lam = lambda_build([A1 A2 A3], r, n);
%   I = eye(2);
%   Lam = [
%     I    O   O   O   O   O
%     O    I   O   O   O   O
%     O    O   I   O   O   O
%     -A3 -A2 -A1  I   O   O
%     O   -A3 -A2 -A1  I   O
%     O    O  -A3 -A2 -A1  I];
  w = lambda_multiply([A1 A2 A3], x, false(r, n));
  ascertain(almostequal(w, Lam*x));
  %
  % LAMBDA_T_MULTIPLY:
  W = lambda_T_multiply([A1 A2 A3], x, r, n);
  ascertain(almostequal(W, Lam'*x));  
  %
  % MORE MISSING VALUES:
  obs = [1 5 6 9 10 12];
  miss = true(r,n); miss(obs) = false;
  Wo = lambda_multiply([A1 A2 A3], x(obs), miss);
  ascertain(almostequal(Wo, Lam(obs,obs)*x(obs)));
  %
  % AND FOR LAMBDA_T_MULTIPLY:
  Wo = lambda_T_multiply([A1 A2 A3], x(obs), r, n, miss);
  ascertain(almostequal(Wo, Lam(obs,obs)'*x(obs))); 
  obs = [1 3 4 5 8 10 11 12];
  miss = true(r,n); 
  miss(obs) = false;
  Wo = lambda_T_multiply([A1 A2 A3], x(obs), r, n, miss);
  ascertain(almostequal(Wo, Lam(obs,obs)'*x(obs))); 
  %
  % NOW TEST WITH EMPTY A:
  w1 = lambda_multiply([], x, false(r, n));
  ascertain(almostequal(w1, x));
  %
  % DERIVATIVE WITH EMPTY A:
  xd = rand(12,1,8);
  [w2,wd] = lambda_multiply([], x, false(r, n), xd);
  ascertain(almostequal(w1, w2));
  ascertain(almostequal(wd, xd));
  %
  % TEST DERIVATIVE FOR COMPLETE DATA AND ZERO xd
  xd = zeros(12,1,8);
  [w, dd] = lambda_multiply([A1 A2], x, false(r, n), xd);
  dd1 = [
    -3  0 -4  0 -1  0 -2  0
     0 -3  0 -4  0 -1  0 -2
    -5  0 -6  0 -3  0 -4  0
     0 -5  0 -6  0 -3  0 -4
    -7  0 -8  0 -5  0 -6  0
     0 -7  0 -8  0 -5  0 -6
    -9  0 -10 0 -7  0 -8  0
     0 -9  0 -10 0 -7  0 -8];
  dd = squeeze(dd(5:12,1,:));
  ascertain(almostequal(dd1, dd));
  %
  % DERIVATIVE FOR COMPLETE DATA AND NONZERO xd
  A = [A1 A2 A3];
  p = 3;
  nPar1 = 0;
  nPar = r^2*p + nPar1;
  xd = rand(length(x), 1, nPar);
  [w, wd] = lambda_multiply(A, x, false(r, n), xd);
  [Lam, Lamd] = lambda_build(A, r, n, nPar);
  for k=1:nPar
    ascertain(almostequal(wd(:,:,k), Lam*xd(:,:,k) + Lamd(:,:,k)*x));
  end
  %
  % NUMERICAL DERIVATIVE COMPARISON
  theta = [A(:); ones(nPar1,1)];
  [d,g,gnum] = diff_test(@wodfun, theta, xd, r, p, n, false(r,n));
  fprintf_if(~qui, 'Derivative test (1): %.1e\n', d);
  ascertain(d<1e-8);  
  %
  % AND MISSING VALUES WITH DERIVATIVE
  x = x(obs);
  xd = xd(obs, :, :);
  [d,g,gnum] = diff_test(@wodfun, theta, xd, r, p, n, miss);
  fprintf_if(~qui, 'Derivative test (2): %.1e\n', d);
  ascertain(d<1e-8);
  %
  disp('OK');
end

function [f,g] = wodfun(theta, xd, r, p, n, miss)
  N = size(xd,1);
  nPar = size(xd, 3);
  x = reshape(xd, N, nPar)*theta;
  A = reshape(theta(1:r^2*p),r,r*p);
  [w, wd] = lambda_multiply(A, x, miss, xd);
  f = w(:);
  g = reshape(wd, N, nPar);
end
