% GEGENBAUERPOLYNOMIAL
%
% Reference: (1) Handbook of Mathematical Functions, Abramowitz and Stegun, p. 782
%
% Input   k        polynomial order
%         LAMBDA   Gegenbauer exponent
%         x[]      vector of points at which to evaluate G_k^LAMBDA(x)
% Output  
%        r[]        recursive caluclation of G_k^LAMBDA(x)
% Called by:
%   1) grp.m, 2) inverseReprojection.m
% Last modified: October 17, 2007


function r = gegenbauerPolynomial(k,LAMDA,x)              
    if k==0 
        r = 1.0;
    elseif (k==1) 
        r = 2.0*LAMDA*x;
    else  % k>1
      k2 = 1.0;          % k-2 
      k1 = 2.0*LAMDA*x;  % k-1
      for i=2:k
         r = (1.0/i)*( 2.0*(i+LAMDA-1.0)*x.*k1 - (i+2.0*LAMDA-2.0).*k2 );
        k2 = k1;
        k1 = r;
      end % for
    end  % else
