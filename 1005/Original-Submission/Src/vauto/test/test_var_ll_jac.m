%TEST_VAR_LL_JAC  Test Jacobian feature of var_ll
%
%  TEST_VAR_LL_JAC compares gradient of missing value log-likelihood given by
%  var_ll with a numerical gradient for a few testcases, when change of
%  variables is used (further testing is carried out by test_var_ll).

function test_var_ll_jac
  fprintf('TESTING JACOBIAN FEATURE OF VAR_LL...')
  rand('state',3); randn('state',1);
  for r=1:2:3
    for p=1:3
      n = p+4;
      X = rand(r,n);
      m1 = r*r*p;      % A
      m2 = r*(r+1)/2;  % Sig
      n1 = ceil(m1/4); % thA
      n2 = 2;          % thSig
      n3 = r;          % mu
      S1 = gallery('lehmer',r);
      S2 = hilb(r); % Sig will be th(n1+1)·S1+th(n1+2)·S2
      Ad = rand(m1, n1);
      s1 = vech(S1);
      s2 = vech(S2);
      Sigd = [s1(:) s2(:)];
      th = (1 + rand(n1 + n2 + n3, 1))*0.2/p/r/r;
      J = blkdiag(Ad, Sigd);
      [A, B, Sig, mu] = theta2mat(th, p, 0, r, J);
      if r>=3 || p >= 3, I = 2; else I = [-4:0 10]; end
      for i = I
        %fprintf('  r,p,i=%d %d %d\n',r,p,i);
        if i==0  % CHECK VAR_LL WITH NO MISSING VALUES
          d = diff_test(@fun_c, th(1:end-r), X, p, r, J);
          if d >= 1e-8, d, end
          ascertain(d < 1e-8);
        end
        Xm = makemissing(X,i);
        d = diff_test(@fun_m, th, Xm, p, r, J);
        if d >= 1e-8, d, end
        ascertain(d < 1e-8)
      end
    end
  end
  disp('  OK')
end

function [ll,lld] = fun_c(theta, X, p, r, J) % Complete data
  [A, B, Sig, mu] = theta2mat(theta, p, 0, r, J);
  if nargout==1, [ll, ok] = var_ll(X, A, Sig);
  else      [ll, ok, lld] = var_ll(X, A, Sig, J); end
  ascertain(ok);
end

function [ll,lld] = fun_m(theta, X, p, r, J) % Missing values
  q = 0;
  miss = isnan(X);
  [A, B, Sig, mu] = theta2mat(theta, p, q, r, J);
  if nargout==1, [ll, ok] = var_ll(X, A, Sig, mu, miss);
  else      [ll, ok, lld] = var_ll(X, A, Sig, mu, miss, J); end
  ascertain(ok);
end
