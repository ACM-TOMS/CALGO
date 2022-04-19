% FRP   Freud reprojection for the removal of Gibbs oscillations
%
% References: (1)  J. of Scientific Computing, 	v. 30, no. 3, (2007) p. 409-440
%             (2)  Applied and Computational Harmonic Analysis, v. 20, no. 1, (2006) 3-25
%  Inputs
%    S     discontinuities (with S[0]=ltBoundaryPt and S[nd+1]=rtBoundaryPt
%   ak     spectral coefficients or function handle
%               spectral coefficients => exact expansion coefficients are calculated
%                     function handle => approximate expansion coefficients are calculated
%   xr     grid to evaluate the FRP approximation on
%  pCh     0 - Fourier, 1 - Chebyshev 
%    N     number of equally spaced function values that are known
%    m     (optional) vector of Freud Polynomial orders for each subinterval         
%  Outputs
%    ug    reprojected approximation
%     a    Freud coefficients in each subinterval
%  Example usage:
%          see fourierFRP_example.m
%  Functions called:
%          1) freudPolynomials.m, 2) offGridFourierInterpolation.m, 3) evaluateChebyshevInterpolant.m
%  Called by:
%          1) postProcessDriver.m, 2) fourierFRP_example.m
%  Notes:
% xi in [-1,1]       x[xi]: [-1,1]->[a,b]
%  x in [a,b]        xi[x]: [a,b]->[-1,1]
%
% Last modified: November 27, 2007

function [ug,a,lambda,m] = frp(S,ak,xr,pCh,N,m)     

    if ~isa(ak,'function_handle'), N = length(ak); end  %  number of spectral coefficients
     M = length(xr);          %  evaluate the GRP approx at M points
    sN = length(S);           %  number of discontinuities and endpoints
   siN = sN - 1;              %  number of sub-intervals
 
  if nargin<6
    for i = 1:siN
      m(i) = round( 0.125*N*( abs(S(i+1)-S(i) )) );                    % (3.17), p. 14, ref (2)
    end
  end
    
% ------------------- coefficient threshold tolerances ---------------------

   if pCh==1             
      tol = 10^(-8);      % Chebshev
   else
      tol = 10^(-12);     % Fourier
   end

% ----------------------------------------------------------------------------
   
    small = 10e-24;
    Nq = N + 50;
    xj = -1 + 2*(0:Nq-1)/Nq;                             %  quadrature (Trapezoid rule) points
    
       mMax = max(m);
       mLimit = 1;
 
       a = zeros(siN,mMax+1);        %  Freud coefficients
       g = zeros(siN,mMax+1);
      uC = zeros(1,Nq);                 %  spectral partial sums
      ug = zeros(1,M);                 %  reconstructed function at the reconstruction pts
 
     if (pCh == 0 & ~isa(ak,'function_handle')),  ak = ((-1).^(0:N-1)).*fftshift(ak)/N; end  % reorder Fourier coefficients
      
 % ------- find Freud coefficients ----------------------------------  
 
     for i = 1:siN                                           % u_N( x[xi] )
     
            sL = S(i); sR = S(i+1); 
            lambda(i) = round( sqrt( 0.5*N*abs(sR-sL) ) - 2*sqrt(2) );            % (3.14), p. 12, ref (2)
            [psi,w,g(i,1:m(i)+1)] = freudPolynomials(xj,m(i),lambda(i),small);
              
            if isa(ak,'function_handle')
            
                 for el=0:m(i)
                     a(i,el+1) = 2*mean( ak(im(sL,sR,xj,0)).*psi(el+1,:).*w );
                  end % el
                  

               for el=2:m(i)       %  adjust m(i) to prevent round-off errors, (4.1)/(4.2), p. 15, ref (2) 
                  if sum(abs(a(i,el-1:el+1)))/3<tol
                     m(i) = max(mLimit,el);
                     break;
                  end
               end  % el
               
               for el=0:m(i)
                 a(i,el+1) = a(i,el+1)./g(i,el+1);
               end % el
            
            
            else
                 if pCh == 0
                   for k=1:Nq
                     uC(k) = offGridFourierInterpolation(ak,im(sL,sR,xj(k),0));                   
                   end
                elseif pCh == 1
                  for k=1:Nq
                    uC(k) = evaluateChebyshevInterpolant(ak,im(sL,sR,xj(k),0));                
                  end
               end
        
        
              for el=0:m(i)
                 a(i,el+1) = 2*mean( uC.*psi(el+1,:).*w );
               end % el
               
               
               for el=2:m(i)
                  if sum(abs(a(i,el-1:el+1)))/3<tol
                     m(i) = max(mLimit,el);
                     break;
                  end
               end  % el
               
               for el=0:m(i)
                 a(i,el+1) = a(i,el+1)./g(i,el+1);
               end % el
              
         end   % else 
  
    end   % i


 % ------- reconstruct function ------------------------------------     
      
         for i=1:M
               for j=1:siN
               
                  if xr(i)>=S(j) & xr(i)< S(j+1)
                     sL = S(j); sR = S(j+1); 
                     
                     for el = 0:m(j)  
                         [fp,notUsed] = freudPolynomials(im(sL,sR,xr(i),1),el,lambda(j),small,g(j,:));
                         ug(i) = ug(i) + a(j,el+1)*fp(el+1,:);                                               % gh*G( xi[x] )
                     end
                     
                   %  break  % get out of the j loop, the interval has been found
                 end  % if in interval
                 
              end  % for j
          end  % for i
          
          if pCh == 1         % for Chebyshev case (xr(M) = 1)/not necessary for Fourier as xr(M)<1
                                        
               ug(M)=0;
               j = siN;
               i = M;
               sL = S(siN); sR = S(siN+1); 
               for el = 0:m(j)  
                  [fp,notUsed] = freudPolynomials(im(sL,sR,xr(i),1),el,lambda(j),small,g(j,:));
                  ug(i) = ug(i) + a(j,el+1)*fp(el+1,:);                               % gh*G( xi[x] )
               end
          
          end
