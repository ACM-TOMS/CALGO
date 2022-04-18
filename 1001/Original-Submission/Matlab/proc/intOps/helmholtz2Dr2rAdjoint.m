%% helmholtz2Dr2rAdjoint
% This function computes the adjoint $V^\ast$ of the volume potential V(f),
%
% $V : L^2(\mathrm{ROI}) \to L^2(\mathrm{ROI})$.
%
% Note that the factor $k^2$ is added in intOpsFuncs and V is redefined.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function uROI = helmholtz2Dr2rAdjoint(fROI, seti)
uROI = restrictCDtoROI(ifft2(seti.kHat'.*fft2(extendROItoCD(fROI,seti.ROImask))),seti.ROImask);
end
