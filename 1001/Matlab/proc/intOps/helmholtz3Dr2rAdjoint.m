%% helmholtz3Dr2rAdjoint
% See the corresponding file for 2D case: <helmholtz2Dr2rAdjoint.html>.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function uROI = helmholtz3Dr2rAdjoint(fROI, seti)
uROI = restrictCDtoROI(ifftn(conj(reshape(seti.kHat,seti.nCD,seti.nCD, seti.nCD)).*fftn(extendROItoCD(fROI,seti.ROImask))),seti.ROImask);
end
