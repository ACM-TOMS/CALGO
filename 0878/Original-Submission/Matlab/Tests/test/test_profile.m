% TEST_PROFILE
%
%   TEST_PROFILE checks the functions profile_ata, profile_mult and
%   profile_back_sub for 9 testcases and 2 missing patterns for each case
%
%   The result of profile_ata(V...) is compared with V'·V and ata_deriv, that of
%   profile_mult(V,Z...) with V'·Z and atb_deriv and profile_back_sub is checked
%   against omega_back_sub. V is calculated as in varma_llm to obtain a
%   convenient test matrix (but could have been assigned random values from
%   scratch).

function test_profile(caseno)
  fprintf('TESTING PROFILE_ATA, PROFILE_MULT AND PROFILE_BACK_SUB... ');
  rand('state',2); randn('state',1);
  if nargin<1, cases=1:9; else cases=caseno; end
  for tcase = cases
    [A, B, Sig, p, q, r, name] = testcase(tcase);
    n = p+q+5;
    nPar = r^2*(p+q) + r*(r+1)/2;
    [C, G, W, S, Cd, Gd, Wd, Sd] = find_CGWS(A, B, Sig);
    [Scol, Scold] = S_extend(A, G, S, n, Gd, Sd);
    [Gneg, Gnegd] = find_Gneg(A, B, C, n, Cd);
    [Su,Olow] = omega_build(S, G, W, p, n);
    [Sud, Olowd] = omega_build_deriv(Sd, Gd, Wd, p, n);
    for k = 0:3
      miss = false(r,n);
      switch k
        case 0, pct = 0;    % Missing percentage
        case 1, pct = 0.25;
        case 2, pct = 0.13;
        case 3, pct = 0.13;
      end
      miss(unique(floor(r*n*rand(round(r*n*pct),1)) + 1)) = true;
      [ko,ro,km] = find_missing_info(miss);
      N = ko(n+1);
      [Suo,Olowo,Suod,Olowod] = omega_remove_miss(Su,Olow,miss,Sud,Olowd);
      [Lu,Ll,info,Lud,Lld] = omega_ltl(Suo, Olowo, p, q, ko, Suod, Olowod);
      [V,Vd] = find_V(G, Gneg, Scol, miss, Gd, Gnegd, Scold);
      Y = profile_back_sub(Lu, Ll, V, ko, km, p, q, q);
      Y1 = omega_back_sub(Lu, Ll, V, p, q, ko);
      ascertain(almostequal(Y,Y1));
      [Y,Yd] = profile_back_sub(Lu, Ll, V, ko, km, p, q, q, Lud, Lld, Vd);
      [Y1,Y1d] = omega_back_sub(Lu, Ll, V, p, q, ko, Lud, Lld, Vd);
      ascertain(almostequal(Y,Y1) && almostequal(Yd,Y1d));
      VV = profile_ata(V, ko, km, p, q);
      ascertain(Lequal(VV, V'*V));
      [VV,VVd] = profile_ata(V, ko, km, p, q, Vd);
      [VV1,VV1d] = ata_deriv(V, Vd);
      ascertain(Lequal(VV,VV1) && Lequal(VVd,VV1d));
      [Acol, Acold] = find_Acol(A, r, nPar);
      [L,Ld] = find_lambda_om(Acol, miss, Acold);
      [mV,nV] = size(V);
      [mL,nL] = size(L);
      V = [V; zeros(N-mV, nV)];
      L = [L; zeros(N-mL, nL)];
      Vd = cat(1, Vd, zeros(N-mV, nV, nPar));
      Ld = cat(1, Ld, zeros(N-mL, nL, nPar));
      Z = rand(N,3);
      Zd = rand(N,3,nPar);
      VZ = profile_mult(V, Z, ko, km, p, q, n);
      LZ = profile_mult(L, Z, ko, km, p, p, n);
      VL = profile_mult(V, L, ko, km, p, q, p);
      LV = profile_mult(L, V  , ko, km, p, p, q);
      ascertain(almostequal(VZ,V'*Z) && almostequal(LZ,L'*Z))
      ascertain(almostequal(VL,V'*L) && almostequal(LV,L'*V))
      [VZ,VZd] = profile_mult(V, Z, ko, km, p, q, n, Vd, Zd);
      [LZ,LZd] = profile_mult(L, Z, ko, km, p, p, n, Ld, Zd);
      [VL,VLd] = profile_mult(V, L, ko, km, p, q, p, Vd, Ld);
      [LV,LVd] = profile_mult(L, V, ko, km, p, p, q, Ld, Vd);
      ascertain(almostequal(VZ,V'*Z) && almostequal(VZd,atb_deriv(V,Vd,Z,Zd)));
      ascertain(almostequal(LZ,L'*Z) && almostequal(LZd,atb_deriv(L,Ld,Z,Zd)));
      ascertain(almostequal(VL,V'*L) && almostequal(VLd,atb_deriv(V,Vd,L,Ld)));
      ascertain(almostequal(LV,L'*V) && almostequal(LVd,atb_deriv(L,Ld,V,Vd)));
    end
  end
  disp('OK');
end

function Lequal = Lequal(x,y)
  Lequal = almostequal(vech(x),vech(y));
end
