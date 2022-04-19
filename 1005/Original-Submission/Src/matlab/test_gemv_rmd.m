% TEST_GEMV_RMD  Check that gemv_rmd works as documented.

function test_gemv_rmd
  for m=1:4
    for n=1:4
      test(m, n);
    end
  end
end

function test(m, n)
  % Test adjoint of mxn matrix times n vector with nf outputs
  alpha = randPM();
  beta = randPM();
  A = randPM(m, n);
  x = randPM(n, 1);
  y = randPM(m, 1);
  u = [vec(A); x; y];
  ftest(@Fmul, @Gmul, u, 'repeats', 2, 'testcase', {m, n, alpha, beta});
end

function v = Fmul(u, m, n, alpha, beta)
  A = reshape(u(1:m*n), m, n);
  x = reshape(u(m*n+1:m*n+n), n, 1);
  y = reshape(u(m*n+n+1:end), m, 1);
  y = alpha*A*x + beta*y;
  v = vec(y);
end

function ua = Gmul(u, ~, ua, va, m, n, alpha, beta)
  nf = size(va,2);
  A = reshape(u(1:m*n), m, n);
  x = reshape(u(m*n+1:m*n+n), n, 1);
  xa = zeros(length(x), nf);
  Aa = zeros([size(A), nf]);
  ya = va;
  [xa,ya,Aa] = gemv_rmd(alpha, A, x, beta, xa, ya, Aa);
  ua = ua + [vec(Aa); xa; ya];
  
  % CHECK THAT UPDATING GIVES THE SAME RESULT AS STARTING WITH ALL ZEROES
  %xi1 = ones(size(x));
  %xi2 = gemvd(a, A, x, b, xi1, ya, Aa);
  %assert(almostequal(xa, xi2));
  
  % Check other calling sequences.
end
