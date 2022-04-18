% SPECTRALFILTER Pth order spectral filters
%
% References: (1) Spectral Methods for Incompressible Viscous Flow, Peyret (2002), ISBN 0-387-95221-7
%             (2) The Erfc-Log Filter and the Asymptotics of the Euler and Vandeven Sum Accelerations", in Proceedings of the
%                 3rd International Conference on Spectral and High Order Methods,  Houston J. Mathematics, Houston, p. 267-276 (1996).
%
% Inputs
%     N,  filter length
%     filterChoice = 1     exponential filter
%                    2     erfc-log
%                    3     Vandeven
%     p,  filter order
%     polyType = 0  Fourier
%                1  Chebyshev
% Outputs
%    sigma  filter coefficients
% Called by:
%   1) filterFourier.m, 2) filterChebyshev.m, 3) postProcessDriver2d.m
% Last modified: October 17, 2007


function sigma = spectralFilter(N,filterChoice,p,polyType)

   switch polyType
     case 1                            %  Chebyshev
          i = 0:N-1;
         nf = N-1;
     otherwise                         %  Fourier
         N2 = floor(0.5*N);
         i(1:N2+1) = 0:N2;
         i(N2+2:2*N2) = -N2+1:-1;
         nf = N2;
   end

   switch filterChoice
        case 1                                                           % exp, b must be even

            sigma = exp( log(eps)*(i/nf).^(p) );                          % reference (1), p. 337, eq. (8.69)

        case 2                                                            % erfc-log
                                                                          % reference (2)
            th = abs(i/nf) - 0.5;
            f1 = 0.5*erfc(0).*(abs(th)<0.0000001);
            warning off
            f2 = (0.5*erfc(2*sqrt(p)*th.*sqrt(-log(1-4*th.^2)./(4*th.^2)))).*(abs(th)>=0.0000001);
            warning on
            f2( find(isnan(f2)) ) = 0;
            sigma = f1 + f2;

        otherwise                                                         % filterChoice=3, Vandeven
                                                                          % reference (1), p. 337, eq. (8.70)
            s = zeros(1,N);
            for j = 0:p-1
               s =  s + ( (-1)^j)*(abs(i/nf).^(2*p-j-1) )/( factorial(p-1-j)*factorial(j)*(2*p-j-1) );
            end
            sigma = 1.0 + ((-1)^p)*(factorial(2*p-1)/factorial(p-1)).*s;

   end
