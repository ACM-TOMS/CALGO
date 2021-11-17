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

disp('EXAMPLE 8: Tries to compute the solution to high accuracy.');
disp('           Although the recursively computed residual is in norm smaller than the tolerance');
disp('           the true residual norm may be larger if a high tolerance is required.');
disp('           Examples 9 and 10 show two remedies.');
disp(' ');
disp('           Sample call:');
disp('           s = 2; tol = 1e-12;');
disp('           [x, flag,relres,iter,resvec] = idrs( A, b, s, tol );');
disp(' ');

tol = 1e-12;

s = 1;
disp('IDR(1) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol );
results;

s = 2;
disp('IDR(2) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol );
results;

s = 4;
disp('IDR(4) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol );
results;

s = 8;
disp('IDR(8) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol );
results;

legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );
hold off;


