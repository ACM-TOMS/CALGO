%% helmholtz3Dc2cf
% See the corresponding file for 2D case: <helmholtz2Dc2cf.html>.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function uCD = helmholtz3Dc2cf(fCD, seti)
uCD = reshape(seti.kHat,seti.nCD,seti.nCD, seti.nCD).*fCD;
end