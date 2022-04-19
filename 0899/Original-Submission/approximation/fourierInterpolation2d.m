%  FOURIERINTERPOLATION2D     2d Fourier interpolation via the FFT
%
% Input:    f[][]             evenly space function values on [-1,1]^2 to be interpolated
%           MX, MY            interpolate to a MX x MY evenly spaced grid on [-1,1]^2
%         
% )utput:   fi[][]            interpolated values 
%           ak[][]            Fourier coefficients
% Called by:
%
% Functions called:
%    1) exampleFunctionSetup2d.m
% Last modified: October 17, 2007


function [fi,ak] = fourierInterpolation2d( f, MX, MY) 

% ---------------------------------------------------------

       [NY,NX] = size( f );                % get input sample size
       gain_x = MX/NX;
       gain_y = MY/NY;

% -----------------------------------------------------------

        f = fft2(f);
        ak = f;
        
% -----------------------------------------------------------
    
    
    if (~mod( NX,2 ) & (MX>NX)) | (mod( NX,2 ) & (MX<NX))
        x_output_space  = max(floor((MX-NX)/2),0) + [1:min(NX,MX)];
        x_input_space   = max(floor((NX-MX)/2),0) + [1:min(NX,MX)];
    else
        x_output_space  = max(ceil((MX-NX)/2),0) + [1:min(NX,MX)];
        x_input_space   = max(ceil((NX-MX)/2),0) + [1:min(NX,MX)];
    end
    
    if (~mod( NY,2 ) & (MY>NY)) | (mod( NY,2 ) & (MY<NY))
       y_output_space  = max(floor((MY-NY)/2),0) + [1:min(NY,MY)];
       y_input_space   = max(floor((NY-MY)/2),0) + [1:min(NY,MY)];
   else
       y_output_space  = max(ceil((MY-NY)/2),0) + [1:min(NY,MY)];
       y_input_space   = max(ceil((NY-MY)/2),0) + [1:min(NY,MY)];
   end
   
% -----------------------------------------------------------
   
    padded_fi    = zeros( MY,MX );                        % perform the up/down sampling
    f            = fftshift(f);
    padded_fi( y_output_space,x_output_space ) = f(y_input_space,x_input_space);
    fi           = (gain_x*gain_y)*ifft2(ifftshift(padded_fi));   
   
% ----------------------------------------------------------- 
  
    fi   = real( fi );
      
    
    
