% FOURIERINTERPOLATION  N equally spaced function values on x = -1 + 2*(0:N-1)/N are
%                       interpolated via the FFT algorithm to M equally spaced points
%                       xp = -1 + 2*(0:M-1)/M
%
%  inputs
%       f   a vector of N evenly spaced function values
%      xp   M evenly spaced points on which to evaluate the interpolant
%  outputs
%     ak    Fourier Coefficients (in Matlab's ordering)
%     fi    interpolated function values on the grid xp
% Last modified: October 17, 2007


function [fi,ak] = fourierInterpolation(f,xp)   

     M = length(xp); 
     N = length(f);    
     gain = M/N; 
     
     f = fft(f);  
     ak = f;        
  
% -----------------------------------------------------------
    
    if (~mod( N,2 ) & (M>N)) | (mod( N,2 ) & (M<N))
        outputSpace  = max(floor((M-N)/2),0) + [1:min(N,M)];
        inputSpace   = max(floor((N-M)/2),0) + [1:min(N,M)];
    else
        outputSpace  = max(ceil((M-N)/2),0) + [1:min(N,M)];
        inputSpace   = max(ceil((N-M)/2),0) + [1:min(N,M)];
    end
    

% -----------------------------------------------------------
   
    padded_fi    = zeros(M,1);                        % perform the up/down sampling
    f            = fftshift(f);
    padded_fi(outputSpace) = f(inputSpace);
    fi           = gain*ifft2(ifftshift(padded_fi));   
   
% ----------------------------------------------------------- 
  
    fi = real( fi );
    fi = fi(:)';
      

    
 
   
   
   
   
