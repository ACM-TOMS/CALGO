%  CHEBYSHEVINTERPOLATION2D 
%
% Inputs
%  u[][]      function values on tensor product CGL grid  
%   xp[]      1d x and y grids from which a 2d grid is formed
%   yp[]         and the interpolant is evaluated on this grid
% Outputs
%   uc[][]    function values at interpolation sites
%   ak[][]    spectral coefficients
% Called by:
%    1) exampleFunctionSetup2d.m
% Functions called:
%    1) evaluateChebyshevInterpolant2d.m,  2) chebyshevCoefficients2d.m, 3) exampleFunctionSetup2d.m
% Last modified: October 17, 2007


  function [uc,ak] = chebyshevInterpolation2d(u,xp,yp)
  
      ak = chebyshevCoefficients2d(u);
      N = length(xp);
      uc = zeros(N,N); 
     
       for k=1:length(xp)
         for j=1:length(yp)
           uc(k,j) = evaluateChebyshevInterpolant2d(ak,xp(k),yp(j));
         end
       end    
  
  
