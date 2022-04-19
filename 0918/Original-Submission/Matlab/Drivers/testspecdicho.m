clear all
close all

n      = 6; 

disp(['     -----------------------'])
disp(['     ---  The matrix  A  ---'])
disp(['     -----------------------'])
A = gallery('frank',n)

pause
disp(['     --------------------------'])
disp(['     ---  Eigenvalues of A  ---'])
disp(['     --------------------------'])
eig(A)

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

specdicho(A)

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

[P,H] = specdicho(A);  
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ]));

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

disp(['     ----------------------------------'])
disp(['     --- Exemple Using the Options  ---'])
disp(['     ----------------------------------'])
opts.mxiter = 20
opts.tol    = 1e-12
opts.r      = 5.38

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause
[P,H] = specdicho(A,opts);
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ]));

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

clear opts
opts.geom      = 'E'

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

[P,H] = specdicho(A,opts); 
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ])); 
disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

opts.geom      = 'I'

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

[P,H] = specdicho(A,opts); 
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ]));  

disp([' '])
disp(['       press a key to continue'])
pause

opts.geom      = 'P'

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

[P,H] = specdicho(A,opts);  
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ])); 

disp(['     ---------------------------'])
disp(['     ---    Pencils Matrix   ---'])
disp(['     ---------------------------'])

B = diag(1:6)

disp(['     ---------------------------------'])
disp(['     --- Unit Circle Using Options ---'])
disp(['     ---------------------------------'])

opts = struct('mxiter',20,'tol',eps,'r',1,'c',0)

disp([' '])
disp(['       press a key to continue'])
pause

[P,H] = specdicho(A,B,opts); 
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ]));  

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

opts.geom      = 'E'

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

[P,H] = specdicho(A,B,opts); 
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ]));  

disp([' '])
disp(['       press a key to continue'])
pause

clear opts
opts.geom      = 'I'

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

[P,H] = specdicho(A,B,opts); 
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ]));  
disp([' '])
disp(['       press a key to continue'])
pause

clear opts
opts.geom      = 'P'

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

[P,H] = specdicho(A,B,opts);  
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ]));

clear all
clc

disp([' '])
disp([' '])
disp(['Example Showing why I might want to use Ho != I'])
disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

A = [1e-3 1e+3;0 1e-3]
B = [1e-5 0;0 1]
[P,H] = specdicho(A,B);
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ]));

disp([' '])
disp(['       press a key to continue'])
disp([' '])
pause

opts.Ho = A'*A + B'*B

[P,H] = specdicho(A,B,opts);
disp(sprintf(['\n NORM(H) = ',num2str(norm(H)) '     NORM(P^2 - P) = ' ...
	      num2str(norm(P^2-P)) '     TRACE(P) = ', num2str(trace(P)) ]));
