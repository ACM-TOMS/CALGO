% FILTERFOURIER2d  Evaluates the spectral filtered 2d fourier approximation
%
% Inputs
%   ak[][]   spectral expansion coefficients
%       M    the filtered approximation is evaluated on a M x M equally spaced tensor product grid on [-1,1]^2
%   sigma    spectral filter
%  Output
%    uf[][]   the filtered approximation
%  Called by:
%          1) postProcessDriver2d.m
% Last modified: October 17, 2007

function uf = filterFourier2d(ak,M,sigma);

% ---------------------------------------------------------

       N = length(sigma);
       gain = M/N;

       akf = ak.*repmat(sigma,N,1).*repmat(sigma',1,N);
% -----------------------------------------------------------


    if (~mod( N,2 ) & (M>N)) | (mod( N,2 ) & (M<N))
        outputSpace  = max(floor((M-N)/2),0) + [1:min(N,M)];
        inputSpace   = max(floor((N-M)/2),0) + [1:min(N,M)];
    else
        outputSpace  = max(ceil((M-N)/2),0) + [1:min(N,M)];
        inputSpace   = max(ceil((N-M)/2),0) + [1:min(N,M)];
    end

% -----------------------------------------------------------

    padded_fi    = zeros( M,M );                       
    f            = fftshift(akf);
    padded_fi( outputSpace,outputSpace ) = f(inputSpace,inputSpace);
    uf           = (gain^2)*ifft2(ifftshift(padded_fi));

% -----------------------------------------------------------

    uf   = real( uf );
