function test_lyap
  % Note: This function currently uses dlyap_octave. It could also
  % be set to use dlyap from Matlab's Control System Toolboxm but 
  % dlyap1 does not work as the test-matrices are not symmetric.
  n = 3;
  A = rand(n);
  M = rand(n);
  x = [vec(A); vec(M)];
  ftest(@Flyap, @Glyap, x, 'tol', 1e-5, 'repeats', 3);
end

function y = Flyap(x)
  global A
  n = sqrt(length(x)/2);
  A = ivec(x(1:n^2));
  M = ivec(x(n^2+1:end));
  X = dlyap(A, M);
  y = vec(X);
end

function xa = Glyap(x, y, xa, ya)
  n = sqrt(length(x)/2);
  A = ivec(x(1:n^2));
  M = ivec(x(n^2+1:end));
  X = ivec(y);
  Xa = ivec(ya);  
  Aa = ivec(xa(1:n^2));
  Ma = ivec(xa(n^2+1:end));
  S = dlyap(A',Xa);
  Ma = Ma + S;
  Aa = Aa + S'*A*X + S*A*X';
  xa = [vec(Aa); vec(Ma)];
end

function A = ivec(u)
  n = sqrt(length(u));
  A = reshape(u, n, n);
end

function S = dlyap(A, X)
  S = dlyap_octave(A, X);
end
