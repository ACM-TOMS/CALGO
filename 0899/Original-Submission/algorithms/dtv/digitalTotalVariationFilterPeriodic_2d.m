% DIGITALTOTALVARIATIONFILTERPERIODIC_2D
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
%        up      the postprocessed function values
% Notes: 
%   uses time-marching (Euler's method) to advance the nonlinear restoration to a steady state 
%   uses a 4 point neighborhood
% Last modified: October 17, 2007

function up = digitalTotalVariationFilterPeriodic_2d(u0,lambda,timeSteps)

     a = 1e-4;
     a = a^2;
     N = length(u0(:,1)) + 1;
     M = length(u0(1,:)) + 1;
     
     U = zeros(N,M);
     U(1:N-1,1:M-1) = u0;
     U(N,:) = U(1,:);                 
     U(:,M) = U(:,1);
     v = U;
     u0 = v;
      
     dt = 0.02;
     
     k = 1;
     s = zeros(N,M);
     
     while k <= timeSteps
     
     i = 2:N-1;
     j = 2:M-1;
     
     s(i,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
     

     i=1;
     s(i,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(N-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
     
     i=N;
     s(i,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(2,j)-U(i,j)).^2 + a );
     
     i = 2:N-1;
     j = 1;
     s(i,j) = sqrt( (U(i,M-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
     
     j = M;
     s(i,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,2)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
     
     i=1; j=1;
     s(i,j) = sqrt( (U(i,M-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(N-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
     
     i=N; j=M;
     s(i,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,2)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(2,j)-U(i,j)).^2 + a );
     
     i=1; j=M;
     s(i,j) = sqrt( (U(i,j-1)-U(i,j)).^2 + (U(i,2)-U(i,j)).^2 + (U(N-1,j)-U(i,j)).^2 + (U(i+1,j)-U(i,j)).^2 + a );
     
     i=N; j=1;
     s(i,j) = sqrt( (U(i,M-1)-U(i,j)).^2 + (U(i,j+1)-U(i,j)).^2 + (U(i-1,j)-U(i,j)).^2 + (U(2,j)-U(i,j)).^2 + a );

  
          i=2:N-1;
          j=2:M-1;
          
            v(i,j) = U(i,j) + dt*( ( U(i+1,j) - U(i,j) ).*( 1 + s(i,j)./s(i+1,j) ) + ...
                                   ( U(i-1,j) - U(i,j) ).*( 1 + s(i,j)./s(i-1,j) ) + ...
                                   ( U(i,j+1) - U(i,j) ).*( 1 + s(i,j)./s(i,j+1) ) + ...
                                   ( U(i,j-1) - U(i,j) ).*( 1 + s(i,j)./s(i,j-1) ) - ...
                                   lambda.*s(i,j).*( U(i,j) - u0(i,j) )                   ); 
                                   
                                  
            i=1; 
            v(i,j) = U(i,j) + dt*( ( U(i+1,j) - U(i,j) ).*( 1 + s(i,j)./s(i+1,j) ) + ...
                                   ( U(N-1,j) - U(i,j) ).*( 1 + s(i,j)./s(N-1,j) ) + ...
                                   ( U(i,j+1) - U(i,j) ).*( 1 + s(i,j)./s(i,j+1) ) + ...
                                   ( U(i,j-1) - U(i,j) ).*( 1 + s(i,j)./s(i,j-1) ) - ...
                                   lambda.*s(i,j).*( U(i,j) - u0(i,j) )                   ); 
                                   
               j=1;
            i=2:N-1;
            v(i,j) = U(i,j) + dt*( ( U(i+1,j) - U(i,j) ).*( 1 + s(i,j)./s(i+1,j) ) + ...
                                   ( U(i-1,j) - U(i,j) ).*( 1 + s(i,j)./s(i-1,j) ) + ...
                                   ( U(i,j+1) - U(i,j) ).*( 1 + s(i,j)./s(i,j+1) ) + ...
                                   ( U(i,M-1) - U(i,j) ).*( 1 + s(i,j)./s(i,M-1) ) - ...
                                   lambda.*s(i,j).*( U(i,j) - u0(i,j) )                   ); 
 
           
            j=M;
            v(i,j) = U(i,j) + dt*( ( U(i+1,j) - U(i,j) ).*( 1 + s(i,j)./s(i+1,j) ) + ...
                                   ( U(N-1,j) - U(i,j) ).*( 1 + s(i,j)./s(N-1,j) ) + ...
                                   ( U(i,2) - U(i,j) ).*( 1 + s(i,j)./s(i,2) ) + ...
                                   ( U(i,j-1) - U(i,j) ).*( 1 + s(i,j)./s(i,j-1) ) - ...
                                   lambda.*s(i,j).*( U(i,j) - u0(i,j) )                   );
                                   
           
   
            v(N,:) = v(1,:);                 
            v(:,M) = v(:,1);



        k = k+1;  U = v; 
   
     end  % while

  up = v(1:N-1,1:M-1);

end



  
  
