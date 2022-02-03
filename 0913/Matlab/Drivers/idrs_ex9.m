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

disp('EXAMPLE 9: Shows the effect of a few replacements of the recursively computed residual by the true residual.');
disp('           Note the difference with example 8.')          
disp(' ');
disp('           s = 2; tol = 1e-12; maxit = 1000;');
disp('           M1 = []; M2 = []; x0 = [];');
disp('           options.replace = 1;');
disp('           [x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );');
disp(' ');

maxit = 1000;
tol = 1e-12;
M1 = [];
M2 = [];
x0 = [];
options.replace = 1;

s = 1;
disp('IDR(1) iteration');
[x, flag,relres,iter,resvec,replacements] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

s = 2;
disp('IDR(2) iteration');
[x, flag,relres,iter,resvec,replacements] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

s = 4;
disp('IDR(4) iteration');
[x, flag,relres,iter,resvec,replacements] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

s = 8;
disp('IDR(8) iteration');
[x, flag,relres,iter,resvec,replacements] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );


