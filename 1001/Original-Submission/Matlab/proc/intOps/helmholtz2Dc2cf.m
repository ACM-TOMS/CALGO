%% helmholtz2Dc2cf
%
% This function computes the volume potential V(f) for the Helmholtz equation 
%
% $V : L^2(\mathrm{CD}) \to L^2(\mathrm{CD})$.
%
% * Note that the factor $k^2$ is added in intOpsFuncs.%
% * Note: function maps column vectors fCD to column vectors uCD.
% * Input and output are Fourier coefficients!
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function uCD = helmholtz2Dc2cf(fCD, seti)

% fROI = fROI(:);
% 'helmholtz2Dr2r'
% size(fROI)
% size(extendROItoCD(fROI,seti.ROImask))
% size(fft2(extendROItoCD(fROI,seti.ROImask)))
% size(reshape(seti.kHat,seti.nCD,seti.nCD))

uCD = reshape(seti.kHat,seti.nCD,seti.nCD).*fCD;
%uROI = uROI(:);

% 'uROI'
% size(uROI)
end
