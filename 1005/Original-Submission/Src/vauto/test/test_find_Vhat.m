% TEST_FIND_VHAT  Check that Vhat is calculated correctly
%
%   TEST_FIND_VHAT checks that the functions find_V and profile_back_sub return
%   the correct Vhat matrix as well as its derivative. Runs testcases 1-7
%
%   TEST_FIND_VHAT(k) tries only testcase k

function test_find_vhat(caseno)
  fprintf('TESTING FIND_V AND PROFILE_BACK_SUB...');
  randn('state',1); rand('state',1);
  if nargin<1, cases=1:7; else cases=caseno; end
  for tcase = cases
    [A, B, Sig, p, q, r, name] = testcase(tcase);
    n = p+q+3;
    [C, G, W, S] = find_CGWS(A, B, Sig);
    Scol = S_extend(A, G, S, n);
    Gneg = find_Gneg(A, B, C, n);
    [Su,Olow] = omega_build(S, G, W, p, n);
    SS = S_build(S, A, G, n);
    Lam = lambda_build(A, r, n);
    %
    % First try with nothing missing:
    ko = 0:r:r*n;
    km = zeros(1,n+1);
    miss = false(r,n);
    [Lu, Ll, info] = omega_factor(Su, Olow, p, q, ko);
    V = find_V(G, Gneg, Scol, miss);
    Vhat = profile_back_sub(Lu, Ll, V, ko, km, p, q, q);
    ascertain(isempty(Vhat));
    %
    % Now some cases with missing values:
    for k = [-4:-1 20 40]
      X = makemissing(ones(r,n), k);
      miss = isnan(X);
      obs = ~miss;
      [ko,ro,km] = find_missing_info(miss);
      [Suo, Olowo] = omega_remove_miss(Su, Olow, miss);
      [Luo, Llo, info] = omega_ltl(Suo, Olowo, p, q, ko);
      V = find_V(G, Gneg, Scol, miss);
      Vhat = profile_back_sub(Luo, Llo, V, ko, km, p, q, q);
      %
      V1 = Lam(obs,obs)*SS(obs,miss) + Lam(obs,miss)*SS(miss,miss);
      Vhat1 = omega_back_sub(Luo, Llo, V1, p, q, ko);
      Vx = zeros(size(Vhat1));
      Vx(1:size(Vhat,1), :) = Vhat;
      ascertain(almostequal(Vx, Vhat1));
      %
      % Finally compare derivatives with numerical ones:
      if k>20 || tcase<4 % Just to save time
        [d,g,gnum] = diff_test(@vhatfun, mat2theta(A, B, Sig), p, q, r, n,miss);
        ascertain(d<1e-8);
      end
    end
  end
  disp('OK');
end

function [f,g] = vhatfun(theta,p,q,r,n,miss);
  [A, B, Sig] = theta2mat(theta,p,q,r);
  nPar = r^2*(p+q) + r*(r+1)/2;
  [ko, ro, km] = find_missing_info(miss); 
  if nargout == 1
    [C,G,W,S] = find_CGWS(A, B, Sig);
    Gneg = find_Gneg(A, B, C, n);
    Scol = S_extend(A, G, S, n);
    [Su,Olow] = omega_build(S, G, W, p, n);
    [Suo, Olowo] = omega_remove_miss(Su, Olow, miss);
    [Luo, Llo, info] = omega_ltl(Suo, Olowo, p, q, ko);
    V = find_V(G, Gneg, Scol, miss);
    Vhat = profile_back_sub(Luo, Llo, V, ko, km, p, q, q);
  else
    [C,G,W,S,Cd,Gd,Wd,Sd] = find_CGWS(A, B, Sig);
    [Scol,Scold] = S_extend(A, G, S, n, Gd, Sd);
    [Gneg, Gnegd] = find_Gneg(A, B, C, n, Cd);
    [Su,Olow] = omega_build(S, G, W, p, n);
    [Sud, Olowd] = omega_build_deriv(Sd, Gd, Wd, p, n);
    [Suo, Olowo,Suod,Olowod] = omega_remove_miss(Su, Olow, miss, Sud,Olowd);
    [Luo, Llo, info,Luod,Llod] = omega_ltl(Suo, Olowo, p, q, ko, Suod, Olowod);
    [V, Vd] = find_V(G, Gneg, Scol, miss, Gd, Gnegd, Scold);
    [Vhat,Vhatd] = profile_back_sub(Luo,Llo, V, ko, km, p, q, q, Luod,Llod, Vd);
    g = reshape(Vhatd,[],nPar);
  end
  f = Vhat(:);
end
