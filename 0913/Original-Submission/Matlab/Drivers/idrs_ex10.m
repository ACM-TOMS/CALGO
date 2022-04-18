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

disp('EXAMPLE 10, If the required tolerance is not met, the process can also be restarted with')
disp('            the found solution as initial guess. This example shows this possibility.')
disp(' ');
disp('           Sample call:');
disp('           s = 2; tol = 1e-12; maxit = 1000;');
disp('           M1 = []; M2 = [];');
disp('           [x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit );');
disp('           if ( flag > 0 ) [x, flag,relres,iter,resvec] = idrs(A,b,s,tol,maxit,M1,M2,x); end;');
disp(' ');

maxit = 1000;
tol = 1e-12;
M1 = [];
M2 = [];

s = 1;
disp('IDR(1) iteration');
[x, flag,relres,iter,resvec] = idrs(A,b,s,tol,maxit);
if ( flag > 0 ) 
   disp(['Flag = ', num2str(flag), '. Required accuracy not achieved, restarting ...'])
   [x, flag, relres, iter_r, resvec_r] = idrs(A,b,s,tol,maxit,M1,M2,x); 
   resvec = [resvec' resvec_r'];
   iter = iter+iter_r;
end
results;

s = 2;
disp('IDR(2) iteration');
[x, flag,relres,iter,resvec] = idrs(A,b,s,tol,maxit);
if ( flag > 0 ) 
   disp(['Flag = ', num2str(flag), '. Required accuracy not achieved, restarting ...'])
   [x, flag,relres,iter_r,resvec_r] = idrs(A,b,s,tol,maxit,M1,M2,x); 
   resvec = [resvec' resvec_r'];
   iter = iter+iter_r;
end
results;

s = 4;
disp('IDR(4) iteration');
[x, flag,relres,iter,resvec] = idrs(A,b,s,tol,maxit);
if ( flag > 0 ) 
   disp(['Flag = ', num2str(flag), '. Required accuracy not achieved, restarting ...'])
   [x, flag,relres,iter_r,resvec_r] = idrs(A,b,s,tol,maxit,M1,M2,x); 
   resvec = [resvec' resvec_r'];
   iter = iter+iter_r;
end
results;

s = 8;
disp('IDR(8) iteration');
[x, flag, relres, iter, resvec] = idrs(A,b,s,tol,maxit);
if ( flag > 0 ) 
   disp(['Flag = ', num2str(flag), '. Required accuracy not achieved, restarting ...'])
   [x, flag,relres,iter_r,resvec_r] = idrs(A,b,s,tol,maxit,M1,M2,x); 
   resvec = [resvec' resvec_r'];
   iter = iter+iter_r;
end
results;

legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );


