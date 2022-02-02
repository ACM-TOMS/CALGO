%  CHEBYSHEVINTERPOLATION
%
%  Inputs
%      u[]   vector of functions values known on the CGL grid on [-1,1]
%     xp[]   grid on [-1,1] to evaluate the interpolant on
%  Outputs
%     ak[]   Chebyshev coefficients
%     uc[]   interpolated function values on the grid xp
%  Functions called:
%   1) chebyshevCoefficients.m, 2) evaluateChebyshevInterpolant.m
%  Called by:
%   1) exampleFunctionSetup1d.m, 2) setUpPdeExample.m
% Last modified: October 17, 2007


   function [uc,ak] = chebyshevInterpolation(u,xp) 
  
        ak = chebyshevCoefficients(u);
     
        for k=1:length(xp)
          uc(k) = evaluateChebyshevInterpolant(ak,xp(k));
        end
  
