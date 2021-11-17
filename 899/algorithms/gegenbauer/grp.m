% GRP   Gegenbauer reprojection for the removal of Gibbs oscillations
%
% References: (1) On the Gibbs phenomenon and its Resolution, SIAM} Review (1997) v. 39, p. 644-668 (and references within)
%
%  Inputs
%    S     discontinuities (with S[0]=ltBoundaryPt and S[nd+1]=rtBoundaryPt
%    L     vector of Gegenbauer Exponents for each subinterval
%    m     vector Gegenbauer Polynomial orders for each subinterval
%   ak     spectral coefficients or function handle
%               spectral coefficients => exact expansion coefficients are calculated
%                     function handle => approximate expansion coefficients are calculated
%   xr     grid to evaluate the GRP approximation on
%  pCh     0 - Fourier, 1 - Chebyshev     
%  Outputs
%    ug    reconstructed data
%    gh    grp coefficients in each subinterval
%  Notes:
% xi in [-1,1]       x[xi]: [-1,1]->[a,b]
%  x in [a,b]        xi[x]: [a,b]->[-1,1]
%  Example usage:
%          see fourierGRP_example.m 
%  Functions called:
%    1) fourierGRP_example.m, 2) offGridFourierInterpolation.m, 3) evaluateChebyshevInterpolant.m
%  Called by:
%          1) postProcessDriver.m, 2) fourierFRP_example.m
% Last modified: October 17, 2007

function [ug,gh] = grp(S,L,m,ak,xr,pCh)     

     if ~isa(ak,'function_handle'), N = length(ak); end          %  number of spectral coefficients
     M = length(xr);          %  evaluate the GRP approx at M points
    sN = length(S);           %  number of discontinuities and endpoints
   siN = sN - 1;              %  number of sub-intervals
 
    Nq = 400;   
    j = 0:Nq;                             %  Chebyshev-Gauss-Lobatto quadrature points
   xj = -cos(j*pi/Nq);                   
 
      gh = zeros(siN,max(m)+1);         %  Gegenbauer coefficients
      uC = zeros(1,Nq+1);               %  Chebyshev partial sums
      ug = zeros(1,M);                  %  reconstructed function at the reconstruction pts

     if (pCh == 0 & ~isa(ak,'function_handle')),  ak = ((-1).^(0:N-1)).*fftshift(ak)/N; end  % reorder Fourier coefficients
      
 % ------- find GRP coefficients ----------------------------------  
 
     wt = pi.*ones(1,Nq+1)./Nq;         % Chebyshev-Gauss-Lobbato quadrature weights
     wt(1) = 0.5*wt(1); 
     wt(end) = 0.5*wt(end); 
     
     for i = 1:siN                                           % u_N( x[xi] )
        sL = S(i); sR = S(i+1); mi = m(i); Li = L(i); 
        
            if isa(ak,'function_handle')
          
                 for el=0:mi
                    hi = ( gamma(Li)*(el+Li))/( sqrt(pi)*gegenbauerPolynomial(el,Li,1.0)*gamma(Li+0.5) );
                    fq = ak(im(sL,sR,xj,0)).*((1.0-xj.^2).^Li).*gegenbauerPolynomial(el,Li,xj);                              %  f( x[xi] )*[(1 - xi^2)^L]*G(xi)
                    gh(i,el+1) = hi*sum(fq.*wt);
                 end % el
        
            else
        
                   if pCh == 0
                       for k=1:Nq+1
                           uC(k) = offGridFourierInterpolation(ak,im(sL,sR,xj(k),0));                 
                       end
                  elseif pCh == 1
                      for k=1:Nq+1
                        uC(k) = evaluateChebyshevInterpolant(ak,im(sL,sR,xj(k),0));                
                      end
                  end
                  
                  for el=0:mi
                      hi = ( gamma(Li)*(el+Li))/( sqrt(pi)*gegenbauerPolynomial(el,Li,1.0)*gamma(Li+0.5) );
                      fq = uC.*((1.0-xj.^2).^Li).*gegenbauerPolynomial(el,Li,xj);                              %  u_N( x[xi] )*[(1 - xi^2)^L]*G(xi)
                      gh(i,el+1) = hi*sum(fq.*wt);
                  end % el
        
          end
    end   % i
       
 % ------- reconstruct function ------------------------------------     
      
         for i=1:M
               for j=1:siN
               
                  if xr(i)>=S(j) & xr(i)< S(j+1)
                     sL = S(j); sR = S(j+1); mi = m(j); Li = L(j);
                     
 
                     for el = 0:mi
                        ug(i) = ug(i) + gh(j,el+1)*gegenbauerPolynomial(el,Li,im(sL,sR,xr(i),1));        % gh*G( xi[x] )
                     end
                     
                     break  % get out of the j loop, the interval has been found
                 end  % if in interval
                 
              end  % for j
          end  % for i
          
 if pCh == 1
          i=M;                                             % for Chebyshev case (xr(M) = 1)/not necessary for Fourier as xr(M)<1
          ug(M)=0;
          sL = S(siN); sR = S(siN+1); mi = m(siN); Li = L(siN);
          for el = 0:mi
            ug(i) = ug(i) + gh(j,el+1)*gegenbauerPolynomial(el,Li,im(sL,sR,xr(i),1));        % gh*G( xi[x] )
          end
 end
