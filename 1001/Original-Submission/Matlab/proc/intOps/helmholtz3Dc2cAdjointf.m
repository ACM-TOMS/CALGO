%% helmholtz3Dc2cAdjointf
% See the corresponding file for 2D case: <helmholtz2Dc2cAdjointf.html>.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function u = helmholtz3Dc2cAdjointf(f, seti)
u = conj(reshape(seti.kHat,seti.nCD,seti.nCD, seti.nCD)).*f;
end
