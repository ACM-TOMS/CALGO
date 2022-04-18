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

disp('EXAMPLE 11, Uses parameters for IDR(s) so that the method is equivalent to Bi-CGSTAB.')
disp(' ');
disp('           s = 1; tol = 1e-8; maxit = 40;');
disp('           M1 = []; M2 = []; x0 = [];');
disp('           options.P = b; options.omega = 0;');
disp('           [x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );');
disp(' ');

maxit = 40;
tol = 1e-8;
M1 = [];
M2 = [];
x0 = [];

s = 1;
options.P = b;
options.omega = 0;
disp('IDR(1) iteration');
[x, flag,relres,iter,resvec] = idrs( A, b, s, tol, maxit, M1, M2, x0, options );
results;

if ( exist('bicgstab') == 2 )
   disp('Bi-CGSTAB iteration');
   [x,flag,relres,iter,resvec_bicgstab] = bicgstab(A,b,tol,maxit/2);
   if ( length(resvec_bicgstab) ==  1+maxit/2 )
      resvec = zeros(maxit+1);
%        bicgstab gives residual norm per whole iteration, not per half! (for use in octave)
      resvec(1:2:maxit+1) = resvec_bicgstab;
      resvec(2:2:maxit) = resvec_bicgstab(1:maxit/2); 
   else
      resvec = resvec_bicgstab;
   end
   results;
   legend('IDR(1)', 'Bi-CGSTAB');
else
   disp('Bi-CGSTAB not implemented!');
   legend('IDR(1)');
end

