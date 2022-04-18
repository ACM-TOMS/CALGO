%% helmholtz3Dr2data
% See the corresponding file for 2D case: <helmholtz2Dr2data.html>.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function measData = helmholtz3Dr2data(fROI, seti)
measData = (seti.measKer*fROI)*seti.dV;
end
