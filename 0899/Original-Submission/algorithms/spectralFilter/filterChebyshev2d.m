% FILTERCHEBYSHEV2d  Evaluates the spectral filtered 2d Chebyshev approximation
%
% Inputs
%    ak   spectral expansion coefficients
%     M   the filtered approximation is evaluated on a M x M tensor product grid on [-1,1]^2
%         defined by the 1d grid xp on [-1,1]
%  sigma  spectral filter
%  Output
%    uf   the filtered approximation
%  Called by:
%          1) postProcessDriver2d.m
%  Functions called:
%          1) evaluateChebyshevInterpolant2d.m
% Last modified: October 17, 2007
%     


  function uf = filterChebyshev2d(ak,xp,sigma);

      N = length(sigma);
      M = length(xp);
      uf = zeros(M,M);


       akf = ak.*repmat(sigma,N,1).*repmat(sigma',1,N);

       for k=1:M
         for j=1:M
           uf(k,j) = evaluateChebyshevInterpolant2d(akf,xp(k),xp(j));
         end
       end
