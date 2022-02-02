%
% IDRS testproblem, to be used with driver test_idrs.m
% Uses results.m, idrs.m, m1.m and m2.m
%
%   Martin van Gijzen
%   Version August 31, 2010
%
%   This software is distributed under the
%   ACM Software Copyright and License Agreement.
%

disp('EXAMPLE 4: SSOR-preconditioned matrix-vector multiplication.')
disp('           Action of preconditioner is implemented using a function.')
disp('           Gives same results as example3.')
disp(' ');
disp('           Sample call:');
disp('           s = 2; tol = []; maxit = [];');
disp('           M1.name = ''m1''; M1.D = diag(diag(A)); M1.L = tril(A)/M1.D; M2.name = ''m2'';  M2.U = triu(A);');
disp('           [x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2 );');
disp(' ');
maxit = [];
tol = [];

M1.name = 'm1';
M1.D = diag(diag(A));
M1.L = tril(A)/M1.D;
M2.name = 'm2';
M2.U = triu(A);

s = 1;
disp('IDR(1) iteration');
[x, flag,relres,iter,resvec] = idrs(A,b,s,tol,maxit,M1,M2 );
results;

s = 2;
disp('IDR(2) iteration');
[x, flag,relres,iter,resvec] = idrs(A,b,s,tol,maxit,M1,M2 );
results;

s = 4;
disp('IDR(4) iteration');
[x, flag,relres,iter,resvec] = idrs(A,b,s,tol,maxit,M1,M2 );
results;

s = 8;
disp('IDR(8) iteration');
[x, flag,relres,iter,resvec] = idrs(A,b,s,tol,maxit,M1,M2 );
results;

legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );


