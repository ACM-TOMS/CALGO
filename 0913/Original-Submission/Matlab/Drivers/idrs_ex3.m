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

disp('EXAMPLE 3: SSOR-preconditioner is passed in the normal way. Preconditioner is decomposed as M1 times M2.'); 
disp('           M1 and M2 are passed as matrices.');
disp(' ');
disp('           Sample call:');
disp('           s = 2; tol = []; maxit = [];');
disp('           [x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2 );');
disp(' ');

tol = [];
maxit = [];

D = diag(diag(A));
M1 = tril(A)/D;
M2 = triu(A);

s = 1;
disp('IDR(1) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2 );
results;

s = 2;
disp('IDR(2) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2 );
results;

s = 4;
disp('IDR(4) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2 );
results;

s = 8;
disp('IDR(8) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2 );
results;

legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );



