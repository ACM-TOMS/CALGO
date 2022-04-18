% OFFGRIDFOURIERINTERPOLATION  Given the Fourier coefficients of a function the DFT is used to
%                              evaluate the interpolant at a single point or on a unequally 
%                              spaced grid.
%
%  Inputs
%        fh[]                   Fourier coefficients
%     x or []                   evaluate the interpolant at x
% N<length(fh)         (optional) use N fh in the appoximation
%  Output
%    s2 or s2[]        the interpolant evaluated at x or x[]
% Called by:
%  1) grp.m, 2) frp.m,  3) inverseReprojection.m
% Last modified: October 17, 2007


   function s2 = offGridFourierInterpolation(fh,x,N)
   
       if nargin<3, N = length(fh); end
   
       s2 = 0;
       for j=1:N
         s2 = s2 + fh(j)*exp(i*(j-N/2-1)*pi*x);
       end
       s2 = real(s2);
 
