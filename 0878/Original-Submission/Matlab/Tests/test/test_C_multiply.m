%TEST_C_MULTIPLY  Tests function C_multiply
%
%  TEST_LAMBDA_MULTIPLY checks that C_multiply returns the right value.

function test_C_multiply
  fprintf('TESTING C_MULTIPLY... ');
  rand('seed',1);
  r = 3;
  n = 4;
  C0 = [1 2 3; 2 3 4; 3 4 5];
  C1 = rand(r);
  C2 = rand(r);
  O = zeros(r);
  C = {C0, C1, C2};
  X = rand(n*r,1);
  miss = false(r,n);
  CC = [C0  C1' C2'  O ;
         O  C0  C1' C2';
         O   O  C0  C1';
         O   O   O  C0];
  Y = C_multiply(C, X, n, miss);
  ascertain(almostequal(Y, CC*X));
  %
  % TRY MISSING VALUES
  miss([1,2,4,8]) = true;
  obs = find(~miss);
  Y = C_multiply(C, X(obs), n, miss);
  ascertain(almostequal(Y, CC(:,obs)*X(obs)));
  fprintf('OK\n');
end
