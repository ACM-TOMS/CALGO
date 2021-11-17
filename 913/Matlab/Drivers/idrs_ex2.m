%
% IDRS testproblem, to be used with driver test_idrs.m
% Uses results.m, idrs.m and mv.m
%
%   Martin van Gijzen
%   Version August 31, 2010
%
%   This software is distributed under the
%   ACM Software Copyright and License Agreement.
%

disp('EXAMPLE 2: SSOR-preconditioned matrix-vector multiplication using Eisenstat trick.')
disp('           The action of the preconditioned MATVEC is implemented using a function.')
disp('           The example also gives comparison with GMRES and Bi-CGSTAB.')
disp('           Note that two-sided preconditioning is applied here. As a result the accuracy')
disp('           in terms of the unpreconditioned residual is poor for all methods.')
disp(' ');
disp('           Sample call:');
disp('           s = 2;');
disp('           K.name = ''mv''; K.D = diag(diag(A)); K.L = tril(A)/K.D; K.U = triu(A)/K.D;');
disp('           f = K.L\b;');
disp('           [x, flag,relres,iter,resvec] = idrs( K, f, s );');
disp('           x = K.D\(K.U\x);');
disp(' ');


K.name = 'mv';
K.D = diag(diag(A));
K.L = tril(A)/K.D;
K.U = triu(A)/K.D;
f = K.L\b;

s = 1;
disp('IDR(1) iteration');
[x, flag,relres,iter,resvec] = idrs( K, f, s );
x = K.D\(K.U\x);
unprec_res = norm(b-A*x)/norm(b);
disp(['Relative unpreconditioned residual norm: ',num2str(unprec_res)]);
results;

s = 2;
disp('IDR(2) iteration');
[x, flag,relres,iter,resvec] = idrs( K, f, s );
x = K.D\(K.U\x);
unprec_res = norm(b-A*x)/norm(b);
disp(['Relative unpreconditioned residual norm: ',num2str(unprec_res)]);
results;

s = 4;
disp('IDR(4) iteration');
[x, flag,relres,iter,resvec] = idrs( K, f, s );
x = K.D\(K.U\x);
unprec_res = norm(b-A*x)/norm(b);
disp(['Relative unpreconditioned residual norm: ',num2str(unprec_res)]);
results;

s = 8;
disp('IDR(8) iteration');
[x, flag,relres,iter,resvec] = idrs( K, f, s );
x = K.D\(K.U\x);
unprec_res = norm(b-A*x)/norm(b);
disp(['Relative unpreconditioned residual norm: ',num2str(unprec_res)]);
results;


if ( exist('bicgstab') == 2  && exist('gmres') == 2 ) 
   tol = 1e-8;
   maxit = min(length(f),1000);
   disp('Bi-CGSTAB iteration');
   [x,flag,relres,iter,resvec] = bicgstab(@mv,f,tol,maxit,[],[],[],K);
   x = K.D\(K.U\x);
   unprec_res = norm(b-A*x)/norm(b);
   disp(['Relative unpreconditioned residual norm: ',num2str(unprec_res)]);
   results;

   disp('GMRES iteration');
   [x,flag,relres,iter,resvec] = gmres(@mv,f,[],tol,maxit,[],[],[],K);
   iter = iter(2);
   x = K.D\(K.U\x);
   unprec_res = norm(b-A*x)/norm(b);
   disp(['Relative unpreconditioned residual norm: ',num2str(unprec_res)]);
   results;
   legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)', 'Bi-CGSTAB', 'GMRES' );
else
   legend('IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );
end




