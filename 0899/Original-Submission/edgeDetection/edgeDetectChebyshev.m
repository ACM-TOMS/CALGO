% EDGEDETECTCHEBYSHEV
%
% References: 1) "Detection of Edges in Spectral Data, Applied and Computational Harmonic Analysis (1999), v. 7, p. 101-135
%             2) "Detection of Edges in Spectral Data II: Nonlinear enhancement", SIAM Journal of Numerical Analysis (2000), v. 38, p. 1389-1408
% Input:
%      ak[]      Chebyshev coefficients
%        J       enhancement threshold
%        Q       enhancement exponent
%       NE      limits how close edges can be together - edges found at
%                 point x(j) will be the strongest discontinuity in the
%                 interval [x(j-NE),x(j+NE)]
%     cfCh      concentration factor choice
%                 0  -  exponential
%                 1  -  linear
% Output:
%        S[]    vector containing the location of the edges
%               with S[1] = -1 and d[sn] = 1.0
% Called by:
%   1) edgeDetectionDriver.m, 2) edgeDectionChebyshev_example.m
% Last modified: July 10, 2008

function [S,uE,uN] = edgeDetectChebyshev(ak,J,Q,NE,cfCh)

   N = length(ak);
   i = 0:N-1;
   x = -cos(pi.*i./(N-1));

% ------------------------------------------------------------------
% --  find edge data, uE -------------------------------------------
% ------------------------------------------------------------------

  invcos = zeros(N,1);
  invcos = acos(x);                % precompute
  const = pi/N;

  uE = zeros(N,1);                 % allocate memorey
  for j=2:N-1;

        k=1:N; k = k';
       
       if cfCh == 1
           sigma = ones(N,1);                                     % linear concentration factor
       else
           sigma = exp(1./(6*k.*(k-1)));  sigma(1)=0;             % exponential concentration factor
       end
       
         temp = ak.*(k-0).*sin( (k-0).*invcos(j) ).*sigma;  
        uE(j) = -sum(temp);

   end

        uE(1) = 0;
        uE(N) = 0;
           uE = const*uE;

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

             for i=1:NE
                if uN(j-i)>0 & uN(j-i)<uN(j)
                  uN(j-i)=0;
                end
                if uN(j+i)>0 & uN(j+i)<uN(j)
                  uN(j+i) = 0;
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
