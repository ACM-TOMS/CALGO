% FILTERFOURIER  spectral filtered fourier approximation
%
% Inputs
%     ak[]   N spectral expansion coefficients
% x or x[]   a single point or vector of points at which the filtered approximation is evaluated
%   fCh      filter choice
%             1     exponential filter
%             2     erfc-log
%             3     Vandeven
%     p   filter order
%  Output
%    fi   the filtered approximation evaluated at x
%  Functions called:
%          1) spectralFilter.m
%  Called by:
%          1) postProcessDriver.m, 2) spectralFilter_example.m
%  Example usage:
%     see spectralFilter_example.m
% Last modified: October 17, 2007
  
  function fi = filterFourier(ak,x,fCh,p);

       sigma = spectralFilter(length(ak),fCh,p,0);

       N = length(ak);
       M = length(x);

       gain = M/N;

       if (~mod( N,2 ) & (M>N)) | (mod( N,2 ) & (M<N))
           outputSpace  = max(floor((M-N)/2),0) + [1:min(N,M)];
           inputSpace   = max(floor((N-M)/2),0) + [1:min(N,M)];
       else
           outputSpace  = max(ceil((M-N)/2),0) + [1:min(N,M)];
           inputSpace   = max(ceil((N-M)/2),0) + [1:min(N,M)];
       end


% -----------------------------------------------------------

    padded_fi    = zeros(M,1);                        % perform the up/down sampling
    ak            = fftshift(sigma.*ak);
    padded_fi(outputSpace) = ak(inputSpace);
    fi           = gain*ifft(ifftshift(padded_fi));

% -----------------------------------------------------------

    fi = real( fi );
    fi = fi(:)';
