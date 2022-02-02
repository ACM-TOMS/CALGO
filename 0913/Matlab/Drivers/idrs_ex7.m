%
% IDRS testproblem, to be used with driver test_idrs.m
% Uses results.m and idrs.m
%
%   Martin van Gijzen
%   Version August 31, 2010
%
%   This software is distributed under the
%   ACM Software Copyright and License Agreement.
%

disp('EXAMPLE 7: Uses complex shadow vectors.');
disp('           Complex shadow vectors speed-up convergence if A is real and strongly nonsymmetric.');
disp('           Effect is strongest for small s, compare with example 1.');
disp('           Of course, complex vectors make the computations more expensive if the problem is real.');
disp(' ');
disp('           Sample call:');
disp('           s = 2; tol = 1e-8; maxit = 1000;');
disp('           M1 = []; M2 = []; x0 = []; im = sqrt(-1);');  
disp('           randn(''state'', 0); P = randn(n,s) + im * randn(n,s); P = orth(P); options.P = P;');
disp('           [x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );');
disp(' ');

maxit = 1000;
tol = 1e-8;
M1 = [];
M2 = [];
x0 = [];
im = sqrt(-1);

s = 1;
randn('state', 0);
P = randn(n,s) + im * randn(n,s);
P = orth(P);
options.P = P;
disp('IDR(1) iteration');
[x, flag, relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

s = 2;
randn('state', 0);
P = randn(n,s) + im * randn(n,s);
P = orth(P);
options.P = P;
disp('IDR(2) iteration');
[x, flag, relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

s = 4;
randn('state', 0);
P = randn(n,s) + im * randn(n,s);
P = orth(P);
options.P = P;
disp('IDR(4) iteration');
[x, flag, relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

s = 8;
randn('state', 0);
P = randn(n,s) + im * randn(n,s);
P = orth(P);
options.P = P;
disp('IDR(8) iteration');
[x, flag, relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );

