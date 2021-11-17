% FILTERCHEBYSHEV  
%
% Inputs
%      ak[]   spectral expansion coefficients
%  x or x[]   a single point or vector of points at which the filtered approximation is evaluated
%   fCh   filter choice
%             1     exponential filter
%             2     erfc-log
%             3     Vandeven
%     p   filter order
%  Output
%    uf   the filtered approximation evaluated at x
%  Functions called:
%          1)  spectralFilter.m,  2) evaluateChebyshevInterpolant.m
%  Called by:
%          1) postProcessDriver.m
% Last modified: October 17, 2007

  function uf = filterChebyshev(ak,x,fCh,p);

       sigma = spectralFilter(length(ak),fCh,p,1);
       akf = ak.*sigma';

       for k=1:length(x)
          uf(k) = evaluateChebyshevInterpolant(akf,x(k));
       end
