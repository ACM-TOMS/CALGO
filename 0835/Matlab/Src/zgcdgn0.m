
function [g,u,v,res] = zgcdgn0(p, q, g0, u0, v0)
%
%  Finds extended GCD of polynomial p and q by Gauss-Newton
%  iteration, such that
%
%         conv(g,u) = p,      conv(g,v) = q
%   
%  Calling syntax:
%    [g,u,v,res] = gcdgn(p, q, g0, u0, v0)
%      
%        INPUT:  p, q -- polynomial coefficients
%                g0, u0, v0 -- initial iterates
%        OUTPUT: g, u, v as described above
%                res
%
   
   m = length(g0)-1; m1 = m+1;    % degree of g (i.e. gcd)
   n = length(u0)-1; n1 = n+1;    % degree of u
   k = length(v0)-1; k1 = k+1;    % degree of v
   lp = length(p);                % length of p
   lq = length(q);                % length of q = p'
   
   % making all polynomials monic
   if p(1) ~= 1, p = p/(p(1)); end;   
   if q(1) ~= 1, q = q/(q(1)); end; 
   if g0(1) ~= 1, g0 = g0/(g0(1)); end; 
   if u0(1) ~= 1, u0 = u0/(u0(1)); end; 
   if v0(1) ~= 1, v0 = v0/(v0(1)); end; 
    
   s = conv(g0,u0)-p; t = conv(g0,v0)-q;
   b = [s(2:lp).'; t(2:lq).'];   
   x = [g0(2:m1).'; u0(2:n1).'; v0(2:k1).'];
   w = ones(length(b),1);
   for j = 2:lp
       if abs(p(j)) > 1 
           w(j-1) = 1/abs(p(j));
       end   
   end;
   
   for j = 2:lq
       if abs(q(j)) > 1 
           w(lp-2+j) = 1/abs(q(j));
       end   
   end;
   
   %bke = [max(abs(b))];
   bke = [max(abs(b.*w))]; 
   
   %fprintf('      g-n %g \n',bke(1));
   
   j = 1;
   while j > 0
       
       A = sylmat(g0,u0,v0);
       
       d = scalsq(A,b,w);
       y = x - d;
       
       g = [1,y(1:m).']; u = [1,y(m+1:m+n).']; v = [1,y(m+n+1:m+n+k).'];
       s = conv(g,u)-p; t = conv(g,v)-q;
       b = [s(2:lp).'; t(2:lq).']; bke = [bke,max(abs(b).*w)]; 
       
       j = j + 1;
       %fprintf('      g-n %g,%g \n',bke(j), norm(d));
       
       if bke(j) >= bke(j-1) 
           g = g0; u = u0; v = v0;
           res = bke(j-1);
           break;
       end;

       g0 = g; u0 = u; v0 = v;
       x = y;
       
   end
