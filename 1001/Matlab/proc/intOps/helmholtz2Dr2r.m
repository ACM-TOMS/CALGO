%% helmholtz2Dr2r
% This function computes the volume potential V(f) for the Helmholtz equation 
% $V : L^2(\mathrm{ROI}) \to L^2(\mathrm{ROI})$
%
% Note: function maps column vectors fROI to column vectors uROI.
%
% Note that the factor $k^2$ is added in intOpsFuncs and V is redefined.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function uROI = helmholtz2Dr2r(fROI, seti)

% fROI = fROI(:);
% 'helmholtz2Dr2r'
% size(fROI)
% size(extendROItoCD(fROI,seti.ROImask))
% size(fft2(extendROItoCD(fROI,seti.ROImask)))
% size(reshape(seti.kHat,seti.nCD,seti.nCD))

uROI = restrictCDtoROI(ifft2(reshape(seti.kHat,seti.nCD,seti.nCD).*fft2(extendROItoCD(fROI,seti.ROImask))),seti.ROImask);
%uROI = uROI(:);

% 'uROI'
% size(uROI)
end
