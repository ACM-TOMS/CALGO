%  EVALUATECHEBYSHEVINTERPOLANT  Clenshaw algorithm to evaluate a Chebyshev interpolant
%                                at each x in O(N) operations
%  Reference:
%   (1) Chebyshev Polynomials, J. Mason and D.Handscomb (2003) ISBN 0849303559, p. 27 
% 
% Input
%   a[]     vector of Chebyshev coefficients
%   x       scalar x to evaluate the interpolant at 
% Output
%   fc      the interpolant evaluated at x
% Called by:
%   1) frp.m, 2) grp.m, 3) filterChebyshev.m
% Last modified: October 17, 2007

  function fc = evaluateChebyshevInterpolant(a,x) 
        
        N = length(a)-1;
        b1 = 0.0;                  
        b0 = 0.0;                  
        b2 = 0.0;
        twox = 2.0*x;          
        for i = N+1:-1:1
            b2 = b1;
            b1 = b0;
            b0 = twox*b1 - b2 + a(i);       %  (2.26)
        end
        
        fc = b0-b1*x;                       %  (2.27)
        
                                 %  Note:
%        fc = 0.5*(b0-b2);       %    (2.28) if a(1) is halved
                                 %     i.e., (sum_1^N)' a_k T_k (x)
       
