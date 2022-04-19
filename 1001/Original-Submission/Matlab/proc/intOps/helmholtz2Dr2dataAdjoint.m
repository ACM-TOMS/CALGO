%% helmholtz2Dr2dataAdjoint
% Computes the adjoint to potential operator that maps functions on ROI to
% data.
%
% If seti.measType = 'farField': evaluates adjoint of far field potential
% $V_\infty$
%
% $V_\infty^\ast : L^2(\bf{S}) \to L^2(\mathrm{ROI})$ 
%  $\quad$ with unit sphere $\bf{S}$.
%
% If seti.measType = 'nearField': evaluates adjoint of volume potential
% from ROI to measPnts (receivers positions):
%
% $V^\ast : L^2(\Gamma_s) \to L^2(\mathrm{ROI})$ $\quad$ with boundary
% $\Gamma_s$.
%
% seti.dS  : approximation of the infinitesimal element on boundary
%
% Note that the factor $k^2$ is added in intOpsFuncs.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
%
function u = helmholtz2Dr2dataAdjoint(f, seti)

u = (seti.measKer')*(f.*seti.dSMeas);

% u = conj(seti.measKer.')*(f.*seti.dSMeas);

% fprintf('\n helmholtz2Dr2dataAdjoint: Max seti.measKer = %5f', max(max(max(abs(seti.measKer)))));
% fprintf('Max ff = %5f', max(max(max(abs(f)))));
% fprintf('seti.dSMeas = %5f \n', max(max(max(seti.dSMeas))));
end
