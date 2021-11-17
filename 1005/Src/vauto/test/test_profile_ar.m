% TEST_PROFILE_AR
%
%   TEST_PROFILE_AR tests the function profile_ar_trisolve with input obtained
%   from omega_ar_factor for pure autoregressive cases from testcase.m and
%   several different missing value patterns for each.
%
%   The results of profile_ar_trisolve is compared with that from
%   omega_ar_trisolve. See also test_profile.
%
%   TEST_PROFILE_AR(N) uses only the N-th testcase.

function test_profile_ar(caseno)
  fprintf('TESTING PROFILE_AR_TRISOLVE... ');
  rand('state',2); randn('state',1);
  if nargin<1, cases=1:11; else cases=caseno; end
  for tcase = cases
    [A, B, Sig, p, q, r, name] = testcase(tcase);
    if q==0
      n = p+5;
      nPar = r^2*p + r*(r+1)/2;
      [C, G, W, S, Cd, Gd, Wd, Sd] = find_CGWS(A, B, Sig);
      [Scol, Scold] = S_extend(A, G, S, n, Gd, Sd);
      [Gneg, Gnegd] = find_Gneg(A, B, C, n, Cd);
      [Su,Olow] = omega_build(S, G, W, p, n);
      [Sud, Olowd] = omega_build_deriv(Sd, Gd, Wd, p, n);
      for k = 0:4
        miss = false(r,n);
        switch k
          case 0, pct = 0;    % Missing percentage
          case 1, pct = 0.25;
          case 2, pct = 0.13;
          case 3, pct = 0.13;
          otherwise pct = 0.20;
        end
        miss(unique(floor(r*n*rand(round(r*n*pct),1)) + 1)) = true;
        [ko,ro,km] = find_missing_info(miss);
        N = ko(n+1);
        [Suo,Olowo,Suod,Olowod] = omega_remove_miss(Su,Olow,miss,Sud,Olowd);
        [V,Vd] = find_V(G, Gneg, Scol, miss, Gd, Gnegd, Scold);
        SigS = mds_set_parmat(p+1, Sig, p+1);
        Sigd = der2array(SigS);
        [Lu,LSig,jpat,Lud,LSigd] = omega_ar_factor(S, Sig, miss, Sd, Sigd{1});
        [Y,Yd] = omega_ar_trisolve(Lu, LSig, V, jpat, ko, 'N', Lud, LSigd, Vd);
        Y1 = profile_ar_trisolve(Lu, LSig, V, jpat, ko, km, p, 0, 'N');
        ascertain(almostequal(Y,Y1));
        [Y1,Y1d]=profile_ar_trisolve(Lu,LSig,V,jpat,ko,km,p,0,'N',Lud,LSigd,Vd);
        ascertain(almostequal(Y,Y1))
        ascertain(almostequal(Yd,Y1d))
      end
    end
  end
  disp('OK');
end
