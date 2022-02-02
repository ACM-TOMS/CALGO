%TEST_OMEGA_FORWARD  Check function omega_forward
%
%  TEST_OMEGA_FORWARD calls omega_forward with a few testcases. Further testing
%  of omega_forward comes with test_omega.
%
function test_omega_forward
  fprintf('PRELIMINARY TEST OF OMEGA_FORWARD AND OMEGA_FORWARD_DERIV...')
  % TEST SIMPLE CASES WITH r=1 AND (p,q) = (1,0), (1,1), (0,1).
  %
  % TEST CASE 1:
  L = [5 0 0; 0 2 0; 0 0 2];
  Lu = 5;
  Ll = [2;2];
  p = 1;
  q = 0;
  nPar = 2;
  nY = 2;
  ko = [0 1 2 3];
  y = rand(3,nY);
  x = omega_forward(Lu, Ll, y, p, q, ko);
  ascertain(almostequal(x, L\y))
  Lud = ones(1,1,nPar);
  Lld = ones(2,1,nPar);
  Lld(:,:,1)=Lld(:,:,1)/2;
  yd = ones(3,nY,nPar);
  xd = omega_forward_deriv(Lu, Ll, Lud, Lld, x, yd, p, q, ko);
  for i=1:nPar
    Ldi = diag([Lud(1,1,i);Lld(:,1,i)]);
    ascertain(almostequal(L*xd(:,:,i), (yd(:,:,i) - Ldi*x)));
  end
  if exist('omega_forward_deriv')==3 % It is a dll-file, so try transpose also
    ydt = permute(yd,[2,1,3]);
    tmin=1; tmax=3;
    xdt = omega_forward_deriv(Lu, Ll, Lud, Lld, x', ydt, p, q,ko,tmin,tmax,'T');
    diff = xd - permute(xdt,[2,1,3]);
    ascertain(norm(diff(:)) < 1e-14);
  end
  % TEST CASE 2:
  q = 1;
  Ll = [1 2;1 2];
  L = [5 0 0; 1 2 0; 0 1 2];
  x = omega_forward(Lu, Ll, y, p, q, ko);
  ascertain(almostequal(x, L\y))
  %
  p = 0;
  x = omega_forward(Lu, Ll, y, p, q, ko);
  ascertain(almostequal(L*x, y))
  %
  % TEST KMIN-KMAX FEATURE:
  y = [1;1];
  x = omega_forward(Lu, Ll, y, p, q, ko, 2);
  ascertain(almostequal([0;x], L\[0;y]));
  x = omega_forward(Lu, Ll, y, p, q, ko, 1, 2);
  ascertain(almostequal(x, L(1:2,1:2)\y))
  y = 3;
  x = omega_forward(Lu, Ll, y, p, q, ko, 1, 1);
  ascertain(almostequal(x, L(1,1)\y));
  x = omega_forward(Lu, Ll, y, p, q, ko, 2, 2);
  ascertain(almostequal(x, L(2,2)\y))
  x = omega_forward(Lu, Ll, y, p, q, ko, 3, 3);
  ascertain(almostequal(x, L(3,3)\y))
  %
  % NOW LET r=2
  r = 2;
  O = zeros(r,r);
  Lu = chol(hilb(r))';
  L1 = ones(r,r);
  L0 = chol(gallery('minij',r))';
  Ll = [L1 L0; L1 L0];
  ko = 0:r:3*r;
  L = [Lu O O; L1 L0 O; O L1 L0];
  nY = 1;
  Y = ones(3*r,nY);
  X = omega_forward(Lu, Ll, Y, p, q, ko); % p=1, q=1
  ascertain(almostequal(X, L\Y))

  % Check derivatives with r=2:
  nPar = 1;
  Lud = tril(ones(size(Lu)));
  L0d = tril(ones(r,r));
  L1d = ones(r,r);
  Lld = [L1d L0d; L1d L0d];
  Ld = [Lud O O; L1d L0d O; O L1d L0d];
  if nPar==2
    Lud = cat(3, Lud, Lud/2);
    Lld = cat(3, Lld, Lld/2);
    Ld = cat(3,Ld,Ld/2);
  end
  Yd = ones(3*r,nY,nPar);
  Xd = omega_forward_deriv(Lu, Ll, Lud, Lld, X, Yd, p, q, ko);
  for i=1:nPar
    ascertain(almostequal(L*Xd(:,:,i), (Yd(:,:,i) - Ld(:,:,i)*X)));
  end
  
  disp('OK');
end
