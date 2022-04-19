% CHEBYSHEVCOEFFICIENTS  Chebyshev coefficients via the FFT
%
% Input
%   f[]     vector of function values on the CGL grid on [-1,1]
% Ouput
%   a[]     coefficients of the interpolating Chebyshev expansion
% Called by:
%  1) chebyshevInterpolation.m
% Last modified: October 17, 2007

 function a = chebyshevCoefficients(f) 
 

     N = length(f);
     f = fliplr(f);   f = f(:);          %   x ordered [1 ... -1], make f corresond
     
     F = ifft([f(1:N);f(N-1:-1:2,:)]);
     a = ([F(1); 2*F(2:N-1); F(N)]);


 
 
