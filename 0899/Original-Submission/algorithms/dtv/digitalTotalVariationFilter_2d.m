% DIGITALTOTALVARIATIONFILTER_2D
%
% References: (1) Digitized {PDE} Method for Data Restoration. Ch. 16, p. 751 to 771, 
%                 in Analytic-Computational Methods in Applied Mathematics (2000), G. Anastassiou editor
%             (2) Digital Total Variation Filtering as Postprocessing for Pseudospectral Methods for Conservation Laws,
%                 Numerical Algorithms, vol. 41, p. 17-33, 2006
% Inputs
%       u0[][]     physical space function values (not spectral coefficients)
%      lambda      fitting parameter
%   timeSteps      number of time marching steps
% Output
%        v[][]     the postprocessed function values
% Called by:
%   1) postProcessDriver2d.m
% Notes: 
%   uses time-marching (Euler's method) to advance the nonlinear restoration to a steady state 
%   uses a 4 point neighborhood
% Last modified: October 17, 2007

function v = digitalTotalVariationFilter_2d(u0,lambda,timeSteps)

     a = 1e-4;
     a = a^2;
     N = length(u0(:,1));
     M = length(u0(1,:));
     v = u0; U = u0;
      
     dt = 0.02;
     
     k = 1;
     s = zeros(N,M);
     
     while k <= timeSteps
     
     i = 2:N-1;
     j = 2:M-1;
  
        s(i,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
        i=1;
        s(1,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 +          0           + (U(i+1,j)-U(i,j)).^2 + a );
        i = 2:N-1; j=M;
        s(i,M) = sqrt( (U(i,j-1)-U(i,j)).^2 +          0           + (U(i-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
        i=N; j = 2:M-1;
        s(N,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 +          0               ) ;
        i = 2:N-1; j=1;
        s(i,1) = sqrt(          0           + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
        i=1; j=1;
        s(1,1) = sqrt(          0           + (U(i,j+1)-U(i,j)).^2 +          0           + (U(i+1,j)-U(i,j)).^2 + a );
        i=N; j=M;
        s(N,M) = sqrt( (U(i,j-1)-U(i,j)).^2 +          0           + (U(i-1,j)-U(i,j)).^2 +          0               );
        i=N; j=1;
        s(N,1) = sqrt(          0           + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 +          0               );
        i=1; j=M;
        s(1,M) = sqrt( (U(i,j-1)-U(i,j)).^2 +          0           +          0           + (U(i+1,j)-U(i,j)).^2 + a );
        
%        v = euler(U,k*dt,dt,@F);


          i=2:N-1;
          j=2:M-1;
          
            v(i,j) = U(i,j) + dt*( ( U(i+1,j) - U(i,j) ).*( 1 + s(i,j)./s(i+1,j) ) + ...
                                   ( U(i-1,j) - U(i,j) ).*( 1 + s(i,j)./s(i-1,j) ) + ...
                                   ( U(i,j+1) - U(i,j) ).*( 1 + s(i,j)./s(i,j+1) ) + ...
                                   ( U(i,j-1) - U(i,j) ).*( 1 + s(i,j)./s(i,j-1) ) - ...
                                   lambda.*s(i,j).*( U(i,j) - u0(i,j) )                   );  
%          end   
%        end

        v(1:N,1) = u0(1:N,1);
        v(1:N,M) = u0(1:N,M);
        v(1,1:M) = u0(1,1:M);
        v(N,1:M) = u0(N,1:M);

        k = k+1;  U = v; 
   
     end  % while
   
% ----------- nested functions ---------------------------------

%{
   function fp = F(u,t)

        fp = zeros(N,M);

        for i=2:N-1
          for j=2:M-1
            fp(i,j) = ( u(i+1,j) - u(i,j) ).*( 1 + s(i,j)./s(i+1,j) ) + ...
                      ( u(i-1,j) - u(i,j) ).*( 1 + s(i,j)./s(i-1,j) ) + ...
                      ( u(i,j+1) - u(i,j) ).*( 1 + s(i,j)./s(i,j+1) ) + ...
                      ( u(i,j-1) - u(i,j) ).*( 1 + s(i,j)./s(i,j-1) ) - ...
                      lambda.*s(i,j).*( u(i,j) - u0(i,j) );  
          end   
        end
      
   end
%}
         
% ----------------------------------------------------------

end



  
  
