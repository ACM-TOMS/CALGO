function test_chol_rmd
  % First call with same parameters as the Fortran program
  % test_spotrf_rmd uses first (to facilitate debugging):
  A = hilb(3);
  L = chol(A,'lower');
  La = gallery('lehmer', 3);
  %Aa = chol_update_rmd(L,La);
  % Then test with five Hilbert matricis:
  for n = [1,2,3,5,10,12]
    %A = hilb(n) + eye(n);
    A = [1 1;1 5];
    x = vech(A);
    ftest(@Fchol, @Gchol, x);
    ftest(@Fchol, @Gchol_update, x);
    ftest(@Fchol, @Gchol1, x);
    ftest(@Fchol, @Gchol2, x);
    chol_blocked(A);
    ftest(@Fchol, @Gchol_blocked, x);
  end
end

function y = Fchol(x)
  LA = ivech(x);
  A = LA + tril(LA,-1)';
  L = chol(A,'lower');
  y = vech(L);
end

function fx = Gchol(x, y, fx, fy)
  fL = ivech(fy);
  L = ivech(y);
  A = ivech(x);
  fA = zeros(size(A));
  fA = chol_rmd(A, L, fA, fL);
  fx = fx + vech(fA);
end

function fx = Gchol_update(~, y, fx, fy)
  % Test version that updates La with Aa (assumes Aa=0, needs only L and La)
  La = ivech(fy);
  L = ivech(y);
  Aa = chol_update_rmd(L,La);
  fx = fx + vech(Aa);
end
function fx = Gchol1(~, y, fx, fy)
  % Test version that updates La with Aa (assumes Aa=0, needs only L and La)
  La = ivech(fy);
  L = ivech(y);
  Aa = chol1_rmd(L,La);
  fx = fx + vech(Aa);
end

function fx = Gchol2(~, y, fx, fy)
  % Test version that updates La with Aa (assumes Aa=0, needs only L and La)
  La = ivech(fy);
  L = ivech(y);
  Aa = chol2_rmd(L,La);
  fx = fx + vech(Aa);
end

function fx = Gchol_blocked(~, y, fx, fy)
  % Test version that updates La with Aa (assumes Aa=0, needs only L and La)
  La = ivech(fy);
  L = ivech(y);
  Aa = chol_blocked_rmd(L,La);
  fx = fx + vech(Aa);
end
