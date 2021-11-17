function test_rotmg
  % Check that rotmg works as documented. The need for this arises because
  % of the bug in the Netlib versions of rotmg (LAPACK v. 3.2-3.8) have.
  x = [2 3]';
  test_one_case([100,0.01], [1/3;-4]);
  test_one_case([4 9], x)
  test_one_case([1 1]*1e-8, x)
  test_one_case([1 1]*1e-32, x)
  test_one_case([1 1]*1e-300, x)
end

function test_one_case(d, x)
  y = [0 0]';
  I = eye(2);
  D = diag(d);
  [dt, y(1), param] = rotmg(d, x);
  Dt = diag(dt);
  H = par2h(param);
  G = Dt^0.5*H/D^0.5;
  assert(almostequal(H*x, y), 'y(2) not 0 when D = diag([4,9])')
  assert(almostequal(G'*G, I))
end
