%% helmholtz3Dr2r
% See the corresponding file for 2D case: <helmholtz2Dr2r.html>.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function uROI = helmholtz3Dr2r(fROI, seti)
uROI = restrictCDtoROI(ifftn(reshape(seti.kHat,seti.nCD,seti.nCD, seti.nCD).*fftn(extendROItoCD(fROI,seti.ROImask))),seti.ROImask);
end