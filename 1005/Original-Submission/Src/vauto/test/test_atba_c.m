% TEST_ATBA_C  Test function atba_c
%
%   TEST_ATBA_C checks that ATBA_C works correctly.

function test_atba_c
  fprintf('TESTING ATBA_C...');
  A1 = rand(4,3);
  A2 = rand(4,3);
  B1 = hilb(4);
  B2 = pascal(4);
  C1 = gallery('lehmer',3);
  C2 = gallery('minij',3);
  x = [0.74; 0.33];
  d = diff_test(@fun, x, A1, A2, B1, B2, C1, C2);
  ascertain(all(d<1e-8));
  A1 = zeros(0,3); A2 = A1;
  d = diff_test(@fun, x, A1, A2, B1, B2, C1, C2);
  %d
  disp('OK');
end

function [f,g] = fun(x,A1,A2,B1,B2,C1,C2)
  [n,m] = size(A1);
  A = x(1)*A1 + x(2)*A2; Ad = cat(3, A1, A2);
  B = x(1)*B1 + x(2)*B2; Bd = cat(3, B1, B2);
  C = x(1)*C1 + x(2)*C2; Cd = cat(3, C1, C2);
  B = B(1:n,1:n);
  Bd = Bd(1:n,1:n,:);
  C = C(1:m,1:m);
  Cd = Cd(1:m,1:m,:);
  D = atba_c(A, B, C);
  [f,h] = atba_c(A, B, C, Ad, Bd, Cd);
  f = vech(f);
  g = vech(h);
  assertEqual(tril(D), tril(A'*B*A + C));
  assertEqual(f, vech(D));
end

function assertEqual(A,B)
  ascertain(isequal(size(A),size(B)) && all(abs(A(:)-B(:)) < 1e-14));
end
