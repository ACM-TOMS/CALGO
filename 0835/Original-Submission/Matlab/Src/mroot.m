function [z,ferr,bkerr,pjcnd,job] = mroot( p, tol, thresh, gamma )
%  MRoot calculates roots and multiplicities of a polynomial p,
%        p(x) = a_1 x^n + a_2 x^n-1 + ... + a_n x + a_n+1,
%  given as an (n+1)-row vector p = (a_1, ..., a_n+1) 
%
%  MRoot is developed by Zhonggang Zeng  (email: zzeng@neiu.edu)
%  This code is freely released for research exchange only. The author 
%  is not responsible for any demage from using this code.
%
%  INPUT :       p -- polynomial coefficients. 
%              tol -- initial residual tolerance (optional); 
%           thresh -- zero singular value shreshold (optional);
%            gamma -- residual growth factor (optional)
%
%  OUTPUT :      z -- roots of polynomial p and multiplicities:
%                      z(:,1) = distinct roots
%                      z(:,2) = corresponding multiplicities
%            pjcnd -- the structure-preserving condition number
%            bkerr -- backward error of z
%             ferr -- the estimated forward error
%              job -- status flag. job = 1  when a
%                     nontrivial multiplicity structure is found,
%                     job = 0 otherwise. If job = 0, the other 
%                     output items should be discarded. 
%            
%  CALL:      >> [z,ferr,berr,cond,job] = multroot(p, tol,thresh,gamma)
%
   % set default parameters if not provided
   if nargin == 1
       tol = 1.e-10; thresh = tol*100; gamma = 100;
   elseif nargin == 2
       thresh = tol*100; gamma = 100;
   elseif nargin == 3
       gamma = 100;
   end;
   %
   % transpose if p is a column vector
   %
   [m,n] = size(p);            
   if m > n, p = p'; end;
   %
   % clear leading/trailing zeros
   %
   n = length(p);
   q = p;
   if q(1) == 0 | q(n) == 0 
       jj = find(p); 
       j1 = min(jj); j2 = max(jj);
       q = p(j1:j2);
       q = q/q(1);
   else
       j1 = 1; j2 = n;
   end;
   %
   % scaling
   %
   m = length(q)-1;
   k = ceil(log(abs(q(m+1)))/(m*log(2)));
   c = 2^(-k); 
   q = q.*(c.^[0:m]);
   %
   [y,bke] = gcdroot(q, tol, thresh, gamma);
   
   if bke < 1.0d-2
       [z,bkerr,pjcnd,job] = pejroot(q,y(:,1).',y(:,2)');

       if job == 1
           if j2 < n
              z = [z;0,n-j2];
           end;
           %
           % show off results
           %
           z(:,1) = z(:,1)/c;
           fprintf('\n');
           fprintf('    !!!THE COMPUTATION IS SUCCESSFUL!!!\n');
           fprintf('\n');
           fprintf('THE PEJORATIVE CONDITION NUMBER:     \t\t\t\t  %g \n',pjcnd);
           fprintf('THE BACKWARD ERROR:                    %6.2e \n',bkerr);
           fprintf('THE ESTIMATED FORWARD ROOT ERROR:      %6.2e \n',...
               2*bkerr*pjcnd);
           fprintf('\n');
           if norm(imag(z(:,1))) == 0 
               fprintf('        computed roots         multiplicities\n');
               fprintf('\n');
               fprintf('%25.15f \t \t \t %3g \n', z');
           else
               fprintf('        computed roots ')
               fprintf('   \t\t\t\t\t\t     multiplicities\n');
               fprintf('\n');
               fprintf('%22.15f + %22.15f i \t \t %3g \n', ...
                   [real(z(:,1)),imag(z(:,1)),z(:,2)]');
           end;
       else
           z = y; bkerr = bke; pjcnd = 1; ferr = 0; job = 0;
       end;
   else
       z = y; bkerr = bke; pjcnd = 1; ferr = 0; job = 0;
   end;
