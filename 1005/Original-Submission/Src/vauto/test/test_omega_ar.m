%TEST_OMEGA_AR  Test omega functions for pure autoregressive case:
%
%  TEST_OMEGA_AR checks the functions omega_ar_factor, omega_ar_trisolve and
%  omega_ar_logdet, both function value and derivative value, by comparing to
%  the values returned by omega_factor, omega_forward and omega_logdet.

function test_omega_ar(varargin)
  [varargin,quiet] = getflags(varargin,'quiet');
  fprintf('TESTING OMEGA_AR_FACTOR, OMEGA_AR_TRISOLVE AND OMEGA_AR_LOGDET...');
  randn('state',1); rand('state',1);
  p = 2; r = 3; n = 6;
  [A,B,Sig,p,q,r]=testcase(p,0,r);
  X = varma_sim(A, B, Sig, n);
  X = makemissing(X, 12);
  miss = isnan(X);
  miss1 = miss(:,1:p);
  miss2 = miss(:,p+1:end);
  PLU = vyw_factorize(A);
  S = vyw_solve(A, PLU, {Sig});
  %
  Su = omega_build(S, {Sig}, {Sig}, p, p); % Test omega_build on finding Su only
  SS = S_build(S, A, {Sig}, p);
  ascertain(almostequal(tril(Su),tril(SS)))
  %
  % Test omega_ar_trisolve, complete data, not transp. (i.e. forward substition)
  ko = 0:r:r*n;
  jpat = ones(1, n-p);
  [Su,Olow] = omega_build(S, {Sig}, {Sig}, p, n);
  [Lu,Ll,info] = omega_factor(Su,Olow,p,q,ko);
  Y = rand(n*r,2);
  Z = omega_forward (Lu, Ll, Y, p, q, ko);
  LSig = {chol(Sig')'};
  Zar = omega_ar_trisolve(Lu, LSig, Y, jpat, ko, 'NoT');
  ascertain(almostequal(Z,Zar));
  %
  % Test omega_ar_trisolve with derivatives, complete data
  nPar = r^2*p + r*(r+1)/2;
  Yd = rand(n*r,2,nPar);
  SigS = mds_set_parmat(p+1, Sig, p+1);
  RHS = vyw_deriv_rhs(A, SigS, S);
  Sd = vyw_solve(A, PLU, RHS);
  Sigd = der2array(SigS);
  LSigd = chol_deriv(LSig{1}, Sigd{1});
  [Sud, Olowd] = omega_build_deriv(Sd, Sigd, Sigd, p, n);
  [Lud, Lld] = omega_factor_deriv(Sud, Olowd, Lu, Ll, p, q, ko);
  Zd = omega_forward_deriv(Lu, Ll, Lud, Lld, Z, Yd, p, q, ko);
  [Zar, Zard] = omega_ar_trisolve(Lu, LSig, Y, jpat, ko, 'NoT', Lud,{LSigd},Yd);
  ascertain(almostequal(Zd,Zard));
  %
  % Check omega_ar_factor Lu:
  [Luo, LSig, jpat] = omega_ar_factor(S, Sig, miss);
  Luo = tril(Luo);
  obs = ~miss1(:);
  Suo = Su(obs,obs);
  ascertain(almostequal(tril(Luo*Luo'), tril(Suo)));
  %
  % Check omega_ar_factor LSig:
  npat = length(LSig);
  for t=p+1:n
    j = jpat(t-p);
    obs = ~miss(:,t);
    ascertain(almostequal(LSig{j}*LSig{j}', Sig(obs,obs)));
  end
  %
  % Check omega_ar_factor derivative:
  theta = mat2theta(A, B, Sig);
  dmax = diff_test(@fun, theta, p, r, miss);
  ascertain(dmax<1e-9);
  %
  % Test omega_ar_trisolve with missing values
  [ko, ro, km] = find_missing_info(miss);
  nobs = ko(n+1);
  Y = rand(nobs,2);
  [Suo,Olowo,Suod,Olowod] = omega_remove_miss (Su, Olow, miss, Sud, Olowd);
  [Luo, Llo, info] = omega_factor (Suo,Olowo,p,q,ko); ascertain(info==0);
  Z = omega_forward (Luo, Llo, Y, p, q, ko);
  Zar = omega_ar_trisolve(Luo, LSig, Y, jpat, ko, 'NoT');
  ascertain(almostequal(Z,Zar));
  %
  % Test omega_ar_trisolve with derivatives
  Yd = rand(nobs,2,nPar);
  [Luod, Llod] = omega_factor_deriv (Suod, Olowod, Luo, Llo, p, q, ko);
  Zd = omega_forward_deriv(Luo, Llo, Luod, Llod, Z, Yd, p, q, ko);
  [Luo,LSig,jpat,Luod,LSigd] = omega_ar_factor(S,Sig,miss,Sd,Sigd{1});
  [Zar,Zard] = omega_ar_trisolve(Luo, LSig, Y, jpat, ko, 'NoT', Luod, LSigd,Yd);
  ascertain(almostequal(Zd,Zard));
  %
  % Check omega_ar_logdet
  [ld, ldd] = omega_logdet(Luo, Llo, p, q, ko, Luod, Llod);
  ldx = omega_ar_logdet(Luo, LSig, jpat);
  [ldar, lddar] = omega_ar_logdet(Luo, LSig, jpat, Luod, LSigd);
  ascertain(almostequal(ld,ldx) && almostequal(ld,ldar) && almostequal(ldd,lddar));
  
  disp(' OK');
end

function [f,g] = fun(theta,p,r,miss)
  nPar = r^2*p + r*(r+1)/2;
  [A, B, Sig] = theta2mat(theta, p, 0, r);
  PLU = vyw_factorize(A);
  S = vyw_solve(A, PLU, {Sig});
  [Lu, LSig, jpat] = omega_ar_factor(S, Sig, miss);
  f = vech(Lu);
  for i=1:length(LSig), f = [f; vech(LSig{i})]; end
  if nargout > 1
    SigS = mds_set_parmat(p+1, Sig, p+1);
    RHS = vyw_deriv_rhs(A, SigS, S);
    Sd = vyw_solve(A, PLU, RHS);
    Sigd = der2array(SigS);
    [Lu,LSig,jpat,Lud,LSigd] = omega_ar_factor(S,Sig,miss,Sd,Sigd{1});
    g = vech(Lud);
    for i=1:length(LSig), g = [g; vech(LSigd{i})]; end
  end
end
