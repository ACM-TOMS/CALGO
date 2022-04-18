function demo_run_groundwater  % Works with Matlab or Octave
  
  % SUMMARY OF THE EXAMPLE MODEL
  % There are gridpoints:
  %
  %   b1----1-----2-----3----b7
  %   |     |x    |  y  |     |
  %   b2----4-----5-----6----b6
  %   |     |     |     |  z  |
  %   +-----b3----b4----b5----+
  %
  % ba are boundary points
  % There is a lake at b1--b2--, 100 m above sea level
  % The sea is at --b6--b7
  % The is a watershed along b1--1--2--3--b7 (no flow across)
  % A river of constant slope runs along --b3--b4--b5--
  % u(i) is the groundwater height at i (i=1...6)
  % There are height measurements at x, y, z of 70, 35, 25 m
  % k(1) is constant x permeability left of 2--5--b4
  % k(2) is constant x permeability right of 2--5--b4
  % There is constant rainfall of 10
  % The resulting Poisson PDE is solved using finite differences 
  % k minimizes squared sum of (modelled-height - measurement)
  
  A{1} = [ % Equation coefficients in left half
    4  -1   0 -2  0  0
    -1  2   0  0 -1  0
    0   0   0  0  0  0
    -1  0   0  4 -1  0
    0 -1/2  0 -1  2  0
    0   0   0  0  0  0];
  A{2} = [ % --in right half
    0   0   0  0  0  0
    0   2  -1  0 -1  0
    0  -1   4  0  0 -2
    0   0   0  0  0  0
    0 -1/2  0  0  2 -1
    0   0  -1  0 -1  4];
  A0 = zeros(6,6);
  h = [100 100 75 50 25 0 0]; % boundary values
  f0 = 10*ones(6,1);          % rain
  f = [h(1) 0 h(7) h(2)+h(3) h(4) h(5)+h(6)];
  F(:,1) = [f(1) f(2)/2  0   f(4) f(5)/2  0  ]';
  F(:,2) = [ 0   f(2)/2 f(3)  0   f(5)/2 f(6)]';
  C = [
    3/8 1/8  0  3/8 1/8   0
    0   1/4 1/4  0  1/4 1/4
    0    0   0   0   0  1/3];
  d = [70, 35, 25 - h(5)/3 - h(6)/3]';
  fun = @(x) groundw_sumsq(x, f0, F, A0, A, C, d);
  opt = optimset('GradObj', 'on', 'Display', 'iter');
  if isempty(ver('Octave')) % Matlab
    opt = optimset(opt, 'LargeScale', 'off', 'DerivativeCheck', 'on');
  end
  k0 = [2;2];
  [theta0, theta0a] = fun(k0);
  fmt1 = 'k%s = (%.3f, %.3f), theta(k%s) = %.3f, k%sa = (%.3f, %.3f)\n';
  [k,fk,~,outp,gk] = fminunc(fun, [2;2], opt); % from Optimization toolbox
  fprintf(fmt1, '0', k0, '0', theta0, '0', theta0a)
  fprintf(fmt1, '*', k, '*', fk, '*', gk)
  fmt2 = 'iterations: %d, function calls: %d\n';
  fprintf(fmt2, outp.iterations, outp.funcCount)
end
