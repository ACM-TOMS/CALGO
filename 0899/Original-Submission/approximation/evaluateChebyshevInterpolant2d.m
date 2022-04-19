% EVALUATECHEBYSHEVINTERPOLANT2D  Clenshaw algorithm to evaluate a Chebyshev interpolant
%                                
%  Reference:
%   (1) Chebyshev Polynomials, J. Mason and D.Handscomb (2003) ISBN 0849303559, p. 27 
% In
%   coeff[][]    Chebyshev coefficients
%   x, y         point (x,y) in [-1,1]^2 to evaluate series at
% Out
%   fc = u_N(x,y) 
% Called by:
%   1) filterChebyshev2d.m, 2) chebyshevInterpolation2d.m
% Last modified: October 17, 2007

  function fc = evaluateChebyshevInterpolant2d(coef,x,y) 
        
   % Clenshaw algorithm to evaluate the sum at each
   % x in O(N) operations
   
    Nx = length(coef(1,:))-1;
    Ny = length(coef(:,1))-1;
    t = zeros(1,Nx);
                                               % ---  evaluate t[m] += aN[m][n]*T(y,n)  via a double summation
    for m=1:Nx
        b1 = 0.0; 
        b0 = 0.0;
        b2 = 0.0; 
        twoy = 2.0*y;
        for i=Ny:-1:1
           b2 = b1; 
           b1 = b0; 
           b0 = twoy*b1 - b2 + coef(m,i); 
        end   
         t(m) = b0-b1*y;                       %  (2.27)
    end
    
    
        b1 = 0.0;              % ---  evaluate u += T(x,m)*t[m]     via a single summation
        b0 = 0.0;  
        b2 = 0.0; 
        twox = 2.0*x;
        for i=Nx:-1:1
           b2 = b1; 
           b1 = b0; 
           b0 = twox*b1 - b2 + t(i); 
        end
        
         fc = b0-b1*x;                       %  (2.27)
