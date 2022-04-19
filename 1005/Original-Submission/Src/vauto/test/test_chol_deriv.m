%TEST_CHOL_DERIV AND FORWARD_SUB_DERIV
%
%  TEST_CHOL_DERIV asserts that the analytical derivatives calculated with
%  chol_deriv and forward_sub_deriv agree with central difference
%  approximations. Use TEST_CHOL_DERIV QUIET to suppress output.
%
function test_chol_deriv(varargin);
  [args,quiet] = getflags(varargin,'quiet');
  fprintf('TESTING CHOL_DERIV AND FORWARD_SUB_DERIV...');
  %
  % FIRST TEST TWO-OUTPUT FEATURES:
  n = 4;
  A = gallery('lehmer',n);
  Ad = cat(3,hilb(n),ones(n));
  [L, Ld] = chol_deriv(A, Ad);
  L = tril(L);
  ascertain(almostequal(L, chol(A)'));
  Ld1 = chol_deriv(L, Ad);
  ascertain(almostequal(Ld, Ld1));
  Ad2 = cat(3,tril(hilb(n)),tril(ones(n)));
  Ld2 = chol_deriv(L,Ad2);
  for i=1:2, Ld(:,:,i)=tril(Ld(:,:,i)); end
  for i=1:2, Ld2(:,:,i)=tril(Ld2(:,:,i)); end
  ascertain(almostequal(Ld, Ld2));
  %
  Y = rand(n,2);
  Yd = rand(n,2,2);
  [X, Xd] = forward_sub_deriv(L, Ld, Y, Yd);
  ascertain(almostequal(L*X, Y));
  Xd1 = forward_sub_deriv(L, Ld, X, Yd);
  ascertain(almostequal(Xd, Xd1));
  %
  % NOW COMPARE ANALYTIC AND NUMERICAL DERIVATIVES:
  theta = [0.2 0.3 0.1];
  d(1) = diff_test(@fun_fact, theta, 1);
  d(2) = diff_test(@fun_fact, theta, n);
  d(3) = diff_test(@fun_forward, theta, 1, 1);
  d(4) = diff_test(@fun_forward, theta, 1, 4);
  d(5) = diff_test(@fun_forward, theta, n, 1);
  d(6) = diff_test(@fun_forward, theta, n, 4);
  ascertain(max(d) < 1e-8)
  disp('  OK')
end
  
function [f,g] = fun_fact(theta, n)
  B1 = gallery('lehmer',n); %...to find some use for these gallery functions
  B2 = repmat(1,n,n);
  B3 = wilkinson(n);
  C = gallery('minij',n);
  p = 1:n;
  A = C + theta(1)*B1 + theta(2)*B2 + theta(3)*B3;
  Ad = cat(3,B1,B2,B3);
  L = chol(A)';
  Ld = chol_deriv(L, Ad);
  f = vech(L);
  g = vech(Ld);
end

function [f,g] = fun_forward(theta, n, m)
  B1 = gallery('lehmer',n); %...to find some use for these gallery functions
  B2 = repmat(1,n,n);
  B3 = wilkinson(n);
  C = gallery('minij',n);
  p = 1:n;
  A = C + theta(1)*B1 + theta(2)*B2 + theta(3)*B3;
  Ad = cat(3,B1,B2,B3);
  L = chol(A)';
  Ld = chol_deriv(L, Ad);
  %
  Y1 = ones(n,m);
  Y2 = repmat(1:m,n,1);
  Y3 = repmat((1:n)',1,m);
  Y = theta(1)*Y1 + theta(2)*Y2 + theta(3)*Y3;
  Yd = cat(3,Y1,Y2,Y3);
  %
  X = L\Y;
  Xd = forward_sub_deriv(L, Ld, X, Yd);  
  f = X(:);
  g = reshape(Xd,n*m,3);
end
