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

disp('EXAMPLE 6: Shows effect of if omega is chosen to minimise the residual norm.');
disp('           The only difference with example 1 is that there omega was determined on basis of accuracy.');
disp(' ');
disp('           Sample call:');
disp('           s = 2; tol = []; maxit = [];');
disp('           M1 = []; M2 = []; x0 = [];');
disp('           options.omega = 0;');
disp('           [x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );');
disp(' ');

options.omega = 0;
maxit = 1000;
tol = 1e-8;
M1 = [];
M2 = [];
x0 = [];

s = 1;
disp('IDR(1) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

s = 2;
disp('IDR(2) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

s = 4;
disp('IDR(4) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

s = 8;
disp('IDR(8) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );
