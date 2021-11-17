% CHEBYSHEVCOEFFICIENTS2d  2d Chebyshev coefficients via the FFT
%
% Input
%   f[][]     vector of function values on the CGL grid on [-1,1]^2
% Ouput
%   B[][]     coefficients of the interpolating Chebyshev expansion on [-1,1]^2
% Called by:
%  1) chebyshevInterpolation2d.m, 2) exampleFunctionSetup2d.m
% Last modified: October 17, 2007

 function B = chebyshevCoefficients2d(f)
 
     f = fliplr(f);            %   x ordered [1 ... -1], make f corresond
     f = flipud(f);            % may need this if function is not symetric
     
    [N,M] = size(f);

    F = ifft([f(1:N,:);f(N-1:-1:2,:)]);
    B = real([F(1,:); 2*F(2:N,:)]);
   
    G = B.';
    F = ifft([G(1:M,:);G(M-1:-1:2,:)]);
    B = real([F(1,:); 2*F(2:M,:)]).';
 
