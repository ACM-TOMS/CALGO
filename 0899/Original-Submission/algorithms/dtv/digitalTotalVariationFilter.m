% DIGITALTOTALVARIATIONFILTER   DTV filtering for non-periodic functions of one variable
%
% References: (1) Digitized {PDE} Method for Data Restoration. Ch. 16, p. 751 to 771, 
%                 in Analytic-Computational Methods in Applied Mathematics (2000), G. Anastassiou editor
%             (2) Digital Total Variation Filtering as Postprocessing for Pseudospectral Methods for Conservation Laws,
%                 Numerical Algorithms, vol. 41, p. 17-33, 2006
% Inputs
%       u0      physical space function values (not spectral coefficients)
%   lambda      fitting parameter
%   timeSteps   number of time marching steps
% Output
%        v      the postprocessed function values
%  Called by:
%          1) postProcessDriver.m
% Notes: 
%   uses time-marching (Euler's method) to advance the nonlinear restoration to a steady state 
% Last modified: October 17, 2007

function v = digitalTotalVariationFilter(u0,lambda,timeSteps)

     a = (1e-4)^2;
     N = length(u0);
     v = u0; U = u0;
      
     dt = 0.02;

     k = 1;
     s = zeros(1,N);
     j = 2:N-1;
     i = 2:N-1;
     while k <= timeSteps
  
        s(1) = sqrt(        0         + (U(2)  -U(1)).^2 + a );
        s(j) = sqrt( (U(j-1)-U(j)).^2 + (U(j+1)-U(j)).^2 + a );
        s(N) = sqrt( (U(N-1)-U(N)).^2 +        0         + a );
  

         v(i) = U(i) + dt*(  ( U(i+1) - U(i) ).*( 1 + s(i)./s(i+1) ) + ...
                             ( U(i-1) - U(i) ).*( 1 + s(i)./s(i-1) ) - ...
                             lambda.*s(i).*( U(i) - u0(i) )               );

        v(N) = u0(N);
        v(1) = u0(1);

        k = k+1;  U = v; 
   
     end  % while
 
end



  
  
