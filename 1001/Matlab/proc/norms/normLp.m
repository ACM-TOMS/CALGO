%% normLp
% Does _not_ compute Schatten norm, but L^p-Norm of measurements
% (equals Schatten norm up to constant if p = 2).
%
%% Syntax
%
%   normA = normLp(A,seti)
%
%% Input Arguments
%
% * seti.pNorm  : p in L^p-Norm
% * seti.dSMeas :   Approximation of the infinitesimal element of a closed contour 
%             with control points.
%             See also <expSetup.html>, <pntsGeometry.html>, and <dS2D.html>.
%
%% Output Arguments
%
% * normA   :   see "More About"
%
%% More About
%
% $\texttt{normA} = \|A\|_P \ \texttt{seti.dSMeas}^{1/P} \quad$
% with $P = \texttt{seti.pNorm}$.
%
%% See Also
% * <normws.html>
% * <addNoise.html> (usage of this function)
%
%% Code
function normA = normLp(A,seti)
normA = norm(A(:),seti.pNorm)*seti.dSMeas^(1/seti.pNorm);
end
