
function [z, bkerr] = gcdroot( p, tol, thresh, growf)  
%  GCDROOT calculates the multiplicity structure and initial root 
%  approximation of a polynomial p given as an (n+1)-row vector 
%                p = (a_1, ..., a_n+1) 
%
%  GcdRoot is developed by Zhonggang Zeng  (email: zzeng@neiu.edu)
%  This code is freely released for research exchange only. The author 
%  is not responsible for any demage from using this code.
%
%  INPUT :       p -- polynomial coefficients. 
%              tol -- initial residual tolerance (optional); 
%           thresh -- zero singular value shreshold (optional);
%            gamma -- residual growth factor (optional)
%
%  OUTPUT :      z -- roots of polynomial p and multiplicities:
%                      z(:,1) = all distinct roots
%                      z(:,2) = corresponding multiplicities
%            bkerr -- backward error of z
%         
%  CALL :   The simplest way to call gcdroot is:   
%               >> z = gcdroot(p)
%           One can also set the 1, 2, or 3 parameters such as
%               >> [z,bkerr] = gcdroot(p,1e-12,1e-8,10);
%               >> z = gcdroot(p,1e-15);
%               >> [z,bke] = gcdroot(p,1e-10,1e-7);
%           etc. 
%           
   
   % set default parameters if not provided
   if nargin == 1
       tol = 1.e-10; thresh = tol*100; gamma = 100;
   elseif nargin == 2
       thresh = tol*100; gamma = 100;
   elseif nargin == 3
       gamma = 100;
   else
       gamma = growf;
   end;
   drop = 5.0d-5;
   
  % make the polynomial monic
   if p(1) == 0, p = p(min(find(p)):n+1); n = length(p)-1; end;
   if p(1) ~= 1, f = p/(p(1)); end; 
   
   n = length(p)-1;    % get polynomial degree 
   q = polyder(p)/n;   % the derevative
   
   ff = p; gg = q;   % back up the polynomial 
   nf = norm(ff,inf);  % the largest coefficient
   
   mx = n; wtol = tol; s0 = 0; s = 0; wthrh = thresh;
   
   k = n;            % the degree of working polynomial
   while k >= 1 
       if k == 1   % the polynomial is linear, GCD = 1
           h = 1; u = f; v = 1; m = 1;
       else
           U = []; R = [];  % initialize input of sylup
           for m = 1:k
               % when m = k, then gcd = 1
               if m == k, h = 1; u = ff; break; end;
               % Update the Sylvester discriminant matrix with QR decomp.
               [U,R] = sylup(ff,gg,m,U,R); s0 = s;
               % Compute the smallest singular value
               [s,x] = uminsv(R,tol);
               if s < wthrh*nf | m == mx | s < drop*s0
                   % a gcd is suspected, get the triplet
                   u0 = x(1:2:2*m+1).'/x(1); v0 = x(2:2:2*m).'/x(2);
                   g0 = lsqdiv(ff,u0);     g0 = g0.'/g0(1);
                   % refine and confirm the gcd
                   [h,u,v,rec] = zgcdgn0(ff, gg, g0, u0, v0);
                   if rec < wtol | m == mx
                       % gcd is confirmed
                       wtol = max([wtol,rec*gamma]); 
                       wthrh = max([wthrh,(s/nf)*100]);
                       break; 
                   end;
               end;
           end;
       end;
       if k == n
           % the first gcd, all distinct roots
           z = roots(u); l = ones(m,1); 
           id = [1:m]; lid = m;
           if m == 1, l = [n];  break; end;
       else
           t = roots(u);
           
           % match roots (avoid matching the same root twice)
           id0 = []; lid0 = 0;
           for j = 1:m
               [s,jj] = min(abs(z(id) - t(j)));
               l(id(jj)) = l(id(jj)) + 1;
               id0 = [id0,id(jj)]; lid0 = lid0+1;
               id(jj) = id(lid); lid = lid-1; id = id(1:lid);
           end;
           id = id0; lid = lid0;
           
           % if there is only one root left, done
           if m == 1, l(id(1)) = l(id(1)) + k - 1; end;
       end;
       if m > 1, k = k - m;  else, k = 0; break;  end;
       if k > 0 
          f = h; g = polyder(f)/k; 
          ff = f; gg = g; nf = norm(ff);  mx = m;
       end;
   end;
   
   ff = [1];  ll = l;    % form the base polynomial
   mm = max(ll); 
   [yy,jj] = sort(abs(z));   y = z(jj); ll = l(jj);
   for kk = 1:mm
       id = find(ll>0);
       ff = conv(ff,poly(y(id)));
       ll = ll - 1;
   end;
   z = [z, l];   % solution

   w = ones(1,length(p));
   for j = 1:length(p), if abs(p(j))>1, w(j)=1/abs(p(j));end,end;
   bkerr = norm( (ff-p).*w , inf ); 
