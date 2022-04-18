%% helmholtz2Dr2data
%
% If seti.measType = 'farField': evaluates far field $V_\infty$ of volume potential V
%
% $V_\infty : L^2(\mathrm{ROI}) \to L^2(\bf{S})$ $\quad$ with unit sphere $\bf{S}$.
%
% If seti.measType = 'nearField': evaluates solution on measPnts (receivers positions)
% via volume potential V:
%
% $V_G : L^2(\mathrm{ROI}) \to L^2(\Gamma_s)$ $\quad$ with boundary
% $\Gamma_s$.
%
% seti.dV : volume of the infinitesimal element on ROI.
%
% Note that the factor $k^2$ is added in intOpsFuncs.
%
%% See Also
%
% * <intOpsFuncs.html>
% * <setMeasKer.html>
%
%% Code
%
function measData = helmholtz2Dr2data(fROI, seti)
measData = (seti.measKer*fROI)*seti.dV;
end
