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

disp('EXAMPLE 1: uses defaults for the input parameters (except for s).'); 
disp('           It also gives comparisons with GMRES and Bi-CGSTAB (if implemented).');
disp(' ');
disp('           Sample call:');
disp('           s = 2;');
disp('           [x, flag,relres,iter,resvec] = idrs( A, b, s );');
disp(' ');

s = 1;
disp('IDR(1) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s );
results;

s = 2;
disp('IDR(2) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s );
results;

s = 4;
disp('IDR(4) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b );
results;

s = 8;
disp('IDR(8) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s );
results;

if ( exist('bicgstab') == 2 && exist('gmres') == 2 )
   tol = 1e-8;
   maxit = min(length(b),1000);

   disp('Bi-CGSTAB iteration');
   [x,flag,relres,iter,resvec] = bicgstab(A,b,tol,maxit);
   results;

   disp('GMRES iteration');
   [x,flag,relres,iter,resvec] = gmres(A,b,[],tol,maxit);
   iter = iter(2);
   results;

   legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)', 'Bi-CGSTAB', 'GMRES' );
else
   legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );
end

