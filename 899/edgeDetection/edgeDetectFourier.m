% EDGEDETECTFOURIER
%
% References: 1) "Detection of Edges in Spectral Data, Applied and Computational Harmonic Analysis (1999), v. 7, p. 101-135
%             2) "Detection of Edges in Spectral Data II: Nonlinear enhancement", SIAM Journal of Numerical Analysis (2000), v. 38, p. 1389-1408
%
% Inputs
%       ak[]    Chebyshev coefficients
%        J      enhancement threshold
%        Q      enhancement exponent
%       NE      limits how close edges can be together - edges found at
%                 point x(j) will be the strongest discontinuity in the
%                 interval [x(j-NE),x(j+NE)]
%     cfCh      concentration factor choice
%                 0  -  exponential
%                 1  -  linear
% Output
%        S[]   vector containing the location of the edges
%              with S[1] = -1 and d[sn] = 1.0
% Called by:
%   1) edgeDetectionDriver.m
% Last modified: July 10, 2008

function [S,uE,uN] = edgeDetectFourier(ak,J,Q,NE,cfCh)

   N = length(ak);
   x = -1 + 2*(0:N-1)/N;   

% ------------------------------------------------------------------
% --  find edge data, uE -------------------------------------------
% ------------------------------------------------------------------

  ak = ak(:)';   % ensure that ak is a row vector                                   
    k=[0:N/2-1 0 -N/2+1:-1];

  if cfCh == 1
     sigma = ones(1,N);                                         % linear concentration factor
  else
     sigma = exp(1./(6*k.*(k-1)));  sigma(isinf(sigma))=0;     % exponential concentration factor
  end 
  
 uE = -sqrt(2)*2*real( ifft(i*pi*k.*ak.*sigma) )/N;
 uE = fliplr(uE);
 
% ------------------------------------------------------------------
% ----- nonlinear enhancement --------------------------------------
% ------------------------------------------------------------------

       uN = zeros(N,1);
       uN = abs(uE);
       j=1:N;
       temp = abs( (N^(0.5*Q))*(uE(j).^Q) );
       I = find( temp < J );
       uN(I) = 0;

% ------------------------------------------------------------------
% ----- check surrounding points -----------------------------------
% ------------------------------------------------------------------

    for j=1:N
       if  (j-NE)>=1 & (j+NE)<=N
           if uN(j)>0                                  % check surronding points

             for k=1:NE
                if uN(j-k)>0 & uN(j-k)<uN(j)
                  uN(j-k)=0;
                end
                if uN(j+k)>0 & uN(j+k)<uN(j)
                  uN(j+k) = 0;
                end
             end % for

         else  % it is within NE grid points of a boundary
            uN(j)=0;
         end
       end % if
  end % j
  
% ------------------------------------------------------------------
% ----- find location of edges -------------------------------------
% ------------------------------------------------------------------

    I = find( uN > eps );
    S(2:2+length(I)-1) = x(I);
    S(1)     =  -1.0;
    S(end+1) =   1.0;

% ------------------------------------------------------------------
