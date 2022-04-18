%TEST_OMEGA Test omega_factor, omega_forward omega_back_sub and omega_logdet
%
%  TEST_OMEGA checks that the functions omega_factor, omega_ltl, omega_forward,
%  omega_back_sub and omega_logdet work correctly for the test-cases that
%  "testcase" provides and no missing values. The missing value case is tested
%  by test_omega_ar and test_varma_llm (and test_var_ll).
%
%  TEST_OMEGA(n) runs only test case # n
%
function test_omega(varargin)
  [varargin, quiet] = getflags(varargin,'quiet');
  fprintf('TESTING OMEGA_FACTOR, OMEGA_FORWARD AND OMEGA_LOGDET...');
  ncase = testcase('number');
  if ~isempty(varargin) cases = varargin{1}; else cases = 1:ncase; end
  for j = cases
    [A, B, Sig, p, q, r, name] = testcase(j);
    fprintf_if(~quiet, '\n  case %-11s r=%d p=%d q=%d ...', name, r, p, q);
    [C, G, W, S] = find_CGWS(A, B, Sig);
    n = p+q+3;
    [Su, Olow] = omega_build(S, G, W, p, n);
    %
    % TEST OMEGA_FACTOR AND OMEGA_LTL (by building a full Omega mat & comparing)
    h = max(p,q);
    [Lu, Ll, info] = omega_factor(Su, Olow, p, q, 0:r:r*n); ascertain(info==0)
    [Lu1, Ll1, info1] = omega_ltl(Su, Olow, p, q, 0:r:r*n); ascertain(info1==0)
    rh = r*h;
    clear LowLeft LowLeftf LowLeftf1 Omega Omegaf Omegaf1
    Omega(1:rh, 1:rh) = Su;
    I = 1:r;
    J = 1:r*(q+1);
    for i=1:n-h
      LowLeft(I,J) = Olow(I,:);
      LowLeftf(I,J) = Ll(I,:);
      LowLeftf1(I,J) = Ll1(I,:);
      I = I + r;
      J = J + r;
    end
    width = min(size(LowLeft,2), r*n);
    Omega(rh+1:r*n,r*n-width+1:r*n) = LowLeft(:,end-width+1:end);
    Omega = tril(Omega) + tril(Omega,-1)';
    [R ,cholinfo] = chol(Omega);
    R1 = flipud(fliplr(chol(flipud(fliplr(Omega)))))';
    if cholinfo==0
      ascertain(cholinfo==0);
      Omegaf (1:rh, 1:rh) = Lu ;
      Omegaf1(1:rh, 1:rh) = Lu1;
      Omegaf (rh+1:r*n,r*n-width+1:r*n) = LowLeftf (:,end-width+1:end);
      Omegaf1(rh+1:r*n,r*n-width+1:r*n) = LowLeftf1(:,end-width+1:end);
      Omegaf  = tril(Omegaf ) + tril(Omegaf ,-1)';
      Omegaf1 = tril(Omegaf1) + tril(Omegaf1,-1)';
      d1 = norm(tril(R' -Omegaf ));
      d2 = norm(tril(R1'-Omegaf1));
      ascertain(d1 < 1e-14);
      ascertain(d2 < 1e-14);
      d = [d1 d2];
      %
      % TEST OMEGA_FORWARD:
      rand('state',0);
      Y = rand(r*n,3);
      X = omega_forward(Lu, Ll, Y, p, q, 0:r:r*n);
      d1 = norm(R'*X - Y);
      ascertain(d1 < 1e-14);
      d(end+1) = d1;
      %
      % TEST OMEGA_BACK_SUB (BOTH EMPTY Y, FULL Y AND ENDING WITH ZEROS):
      X = omega_back_sub(Lu, Ll, zeros(r*n,0), p, q, 0:r:r*n);
      ascertain(isequal(size(X), [r*n,0]));
      X = omega_back_sub(Lu, Ll, zeros(0,3), p, q, 0:r:r*n);
      ascertain(isequal(size(X), [0,3]));
      Y1 = Y(1:r*(n-2), :); Y1z = [Y1; zeros(r*2,3)];
      Y2 = Y(1:r*2, :);     Y2z = [Y2; zeros(r*(n-2),3)];
      X = omega_back_sub(Lu, Ll, Y, p, q, 0:r:r*n);
      X1 = omega_back_sub(Lu, Ll, Y1, p, q, 0:r:r*n);
      X2 = omega_back_sub(Lu, Ll, Y2, p, q, 0:r:r*n);
      d1 = norm(R*X - Y);
      X1z = [X1; zeros(r*2,3)];
      X2z = [X2; zeros(r*(n-2),3)];
      d2 = norm(R*X1z - Y1z);
      d3 = norm(R*X2z - Y2z);
      ascertain(d1 < 1e-14);
      ascertain(d2 < 1e-14);
      ascertain(d3 < 1e-14);
      d = [d d1 d2 d3];
      %
      % TEST OMEGA_FORWARD WHEN Y STARTS WITH ZEROES
      Y = [rand(r*(n-1),3)];
      X = omega_forward(Lu, Ll, Y, p, q, 0:r:r*n, 2);
      O = zeros(r,3);
      d1 = R'*[O; X] - [O; Y];
      ascertain(norm(d1)<1e-14);
      d(end+1) = norm(d1);
      %
      % TEST OMEGA_LOGDET:
      ld = omega_logdet(Lu, Ll, p, q, 0:r:r*n);
      d1 = abs((ld - log(det(Omega)))/ld);
      ascertain(d1 < 5e-15);
      
      fprintf_if(~quiet, '  max diff=%.1e', max(d));
    else
      fprintf_if(~quiet, '  Omega not positive definite');
    end
  end
  if length(cases)==1, disp('  OK'); return, end
  %
  % TEST NON-POSITIVE-DEFINITE FAILURE
  [A, B, Sig, p, q, r, name] = testcase('mediumARMA1');
  n = 12;
  [C, G, W, S] = find_CGWS(A, B, Sig);
  [Su, Olow] = omega_build(S, G, W, p, n);
  Olow(end,end) = -1;
  [Lu, Ll, info] = omega_factor(Su, Olow, p, q, 0:r:r*n);
  ascertain(info == n*r);
  Su(3,3) = -1;
  [Lu, Ll, info] = omega_factor(Su, Olow, p, q, 0:r:r*n);
  ascertain(info == 3);
  fprintf_if(~quiet, '\n');
  disp('  OK');
end
