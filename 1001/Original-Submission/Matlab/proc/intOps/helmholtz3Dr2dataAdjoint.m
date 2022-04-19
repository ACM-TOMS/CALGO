%% helmholtz3Dr2dataAdjoint
% See the corresponding file for 2D case: <helmholtz2Dr2dataAdjoint.html>.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function u = helmholtz3Dr2dataAdjoint(f, seti)
u = (seti.measKer')*(f.*seti.dSMeas);
% u = conj(seti.measKer.')*(f.*seti.dSMeas);
end

