%% helmholtz2Dc2cAdjointf
%
% This function computes the adjoint $V^*$ of the volume potential V(f),
%
% $V : L^2(\mathrm{ROI}) \to L^2(\mathrm{ROI})$,
%
% using Fourier coefficients as input and output.
%
% Note that the factor $k^2$ is added in intOpsFuncs.
%
%% See Also
%
% * <intOpsFuncs.html>
%
%% Code
function u = helmholtz2Dc2cAdjointf(f, seti)
u = conj(reshape(seti.kHat,seti.nCD,seti.nCD)).*f;
end
