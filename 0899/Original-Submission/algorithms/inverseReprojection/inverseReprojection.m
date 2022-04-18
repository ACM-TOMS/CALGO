% INVERSEREPROJECTION  Inverse Polynomial Reprojection for Fourier Approximations
%
% References: (1) J. of Compuational and Applied Mathematics, 161 (2003), 41-65
%             (2) J. of Compuational and Applied Mathematics, 172 (2004), 131-151
%             (3) J. of Compuational and Applied Mathematics, 170 (2004), 303-315
%             (4) J. of Scientific Computing, 25, no. 3, December 2005, 367-399
% Inputs
%   S                 edge locations
%   ML                vector Gegenbauer Polynomial orders for each subinterval
%   xr                grid to evaluate the reconstruction on
%    A                either a vector of Fourier coefficients or else a handle to a function
%   Lambda (optional) Gegenbauer parameter, only affects cond(W)
%   Nq     (optional)  number of quadrature points
% Outputs
%   up                 the reconstructed function
%  condW               condition number of the transformation matrix W
% Functions called:
%    1) gegenbauerPolynomial.m, 2) offGridFourierInterpolation.m
%  Called by:
%          1) postProcessDriver.m
% Notes - programmed for the following cases
%   1) f known at grid points
%   2) approximate Fourier coefficients known
% Last modified: October 17, 2007


function [up,condW] = inverseReprojection(S,ML,xr,A,Lambda,Nq)

     if nargin<6
        Nq=600;
        if nargin<5, Lambda=0.5; end         
     end

    tic
    
      N = sum(ML);            % number of Fourier coefficients used in the reprojection
     nh = (N-1)/2;
    
     sN = length(S);          %  number of discontinuities and endpoints
    siN = sN - 1;             %  number of sub-intervals
    
    if ~isa(A,'function_handle')
       N2 = length(A);
       ak = ((-1).^(0:N2-1)).*fftshift(A)/N2;
    end
    
                          
      [w,xj] = gaussLegendreQuadratureNodesAndWeights(Nq);     %  Gauss Legendre quadrature nodes and weights

% ------------ elements of the transformation matrix via Gauss-Legndre quadrature ----------------------------------
        
         W = zeros(N,N);                                     %  transformation matrix

         mSum = 0;
         for j = 1:siN                                          
            sL = S(j); sR = S(j+1); mi = ML(j); 
            ck = 0;
            for k=-nh:nh  
              ck = ck + 1;
              colNumber = mSum;
              for el=0:mi-1
                     colNumber = colNumber + 1;
                     fq = 0.5*exp(-i.*k.*pi.*im(sL,sR,xj,0)).*gegenbauerPolynomial(el,Lambda,xj);    % (39) ref 4, p. 380
                     W(ck,colNumber) = 0.5*abs(sR-sL)*dot(w,fq); 
              end %l
            end  %k
            mSum = mSum + mi;
        end  %j    
            
            
        condW = cond(W);

    
% ---- calculate the N Fourier Coefficients used in the projection via Gauss Legendre quadrature -------------------
% -- Note: the Fourier coefficients and ellements of W must be calculated "consistently".  See ref. (2), p. 146 ----

    fh = zeros(N,1);                     
    ct = 1;
    
        if isa(A,'function_handle')
            for k=-nh:nh    
                   fq = 0.5*A(xj).*exp(-i*k*pi.*xj);
               fh(ct) = dot(w,fq);  
                   ct = ct + 1;
            end
        else                                     % the underlying function is unknown - only the Fourier coeffients are know
            for k=-nh:nh                         % such as would be the case in a pseudospectral approximation
               fq = 0.5*offGridFourierInterpolation(ak,xj).*exp(-i*k*pi.*xj);
               fh(ct) = dot(w,fq);  
               ct = ct + 1;
            end
        end
    
    
% ---------------------------------------------------------------------------------------
         
      g = W\fh;                %  N Gegenbauer coefficients, (42) ref 4, p. 380
   
   
% ---------- evaluate the reprojection ----------------------- 

    M = length(xr);
    up = zeros(1,M);
    mSum = 0; 
    
    for k=1:M
         
               for j=1:siN
               
                  if xr(k)>=S(j) & xr(k)< S(j+1)
                     sL = S(j); sR = S(j+1); 
                     
                     mi = ML(j);
                     
                     if j==1, mSum=0; else
                     mSum = sum(ML(1:j-1)); end
                     
                     for el = 0:mi-1                                                                           % (36) ref 4, p. 380
                        up(k) = up(k) + g(mSum+el+1)*gegenbauerPolynomial(el,Lambda,im(sL,sR,xr(k),1));        % gh*G( xi[x] )
                     end
                     
                     break  % get out of the j loop, the interval has been found
                 end  % if in interval
                 
              end  % for j
              
    end % k
    
    up = real(up);
   
% --------------- subfunction ------------------------------  

   
    
 function [w,xj] = gaussLegendreQuadratureNodesAndWeights(Nq)
 
        beta = .5./sqrt(1-(2*(1:Nq-1)).^(-2));
        [V,D] = eig(diag(beta,1) + diag(beta,-1));
         xj = diag(D); [xj,i] = sort(xj);
         w = 2*V(1,i).^2;
