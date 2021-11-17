%% intOpsFuncs
% Construct function handles on various integral operators.
%
%% Syntax
%
%   [V, S, VGStar, VStar, TStar, QU, QStarU, Vqf, VqStarf] = intOpsFuncs(qROI, seti)
%
%% Description
% |[V, S, VGStar, VStar, TStar, QU, QStarU, Vqf, VqStarf] = intOpsFuncs(qROI, seti)|
% constructs several function handles on various integral operators 
% (in case of 2D and 3D scattering problem).
%
%% Input Arguments
%
% * qROI :  contrast in ROI
% * seti :  structural array
%
%% Output Arguments
%
% Note that the suffix "Star" is used to mark adjoints.
%
% * V       : volume potential $V: L^2(\mathrm{ROI}) \to L^2(\mathrm{ROI})$,
%             defined in <solveLippmannSchwinger.html> in "More About".
% * S       : S = @(s) simo(qROI, s, seti), see <simo.html>.
% * VGStar  : $L^2(\Gamma_s) \to L^2(\mathrm{ROI})$
%             (receivers are placed on manifold $\Gamma_s$.)
% * VStar   : $V^\ast: L^2(\mathrm{ROI}) \to L^2(\mathrm{ROI})$
% * TStar   : $T(q)^\ast : L^2(\mathrm{ROI}) \to L^2(\mathrm{ROI})$
% * QU      : |QU = @(x) qROI .* x;|
% * QStarU  : |QStarU = @(x) conj(qROI) .* x;|
% * Vqf     : $V: L^2(\mathrm{CD}) \to L^2(\mathrm{CD})$ with respect to
%             contrast q
% * VqStarf : $V: L^2(\mathrm{CD}) \to L^2(\mathrm{CD})$ with respect to
%             complex conjugate contrast q.
%
%% More About
%
% For reference see [1] or [2, Sec. 3].
%
% * QU, QStarU : workaround for helmholtzHMode. 
%  The other models often need |qROI.*x| or |conj(qROI).*x| where HMode only needs |x|.
%  (Note that helmholtzHMode is not supported in public version.)
%
%% References
%
% * [1] David Colton and Rainer Kress. _Inverse Acoustic and Electromagnetic Scattering Theory_. Springer, New York, 2013.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <simo.html>
% * <solveLippmannSchwinger.html>
%
% _Corresponding functions for helmholtz2D:_
%
% * <helmholtz2Dc2cAdjointf.html>
% * <helmholtz2Dc2cf.html>
% * <helmholtz2Dr2data.html>
% * <helmholtz2Dr2dataAdjoint.html>
% * <helmholtz2Dr2r.html>
% * <helmholtz2Dr2rAdjoint.html>
%
% _Corresponding functions for helmholtz3D:_
%
% * <helmholtz3Dc2cAdjointf.html>
% * <helmholtz3Dc2cf.html>
% * <helmholtz3Dr2data.html>
% * <helmholtz3Dr2dataAdjoint.html>
% * <helmholtz3Dr2r.html>
% * <helmholtz3Dr2rAdjoint.html>
%
%% Code
function [V, S, VGStar, VStar, TStar, QU, QStarU, Vqf, VqStarf] = intOpsFuncs(qROI, seti)

if strcmp(seti.model, 'helmholtz2D')
    V  = @(x) seti.k^2.*helmholtz2Dr2r(x, seti); % volume potential V: L^2(ROI) -> L^2(ROI)

    % Vf: volume potential V essentially from  L^2(CD) to L^2(CD), 
    %     but with Fourier coefficients as input and output,
    %     see Vqf for correct usage.
    Vf  = @(x) seti.k^2.*helmholtz2Dc2cf(x, seti);
    
    % Vqf: respects contrast q and that Vf has Fourier coefficients as
    % input and output, 
    Vqf = @(x) Vf(fftn(extendROItoCD(qROI .* restrictCDtoROI(ifftn(x),seti.ROImask),seti.ROImask)));
    
    VGStar = @(x) seti.k^2.*helmholtz2Dr2dataAdjoint(x, seti);
    
    VStar  = @(x) seti.k^2.*helmholtz2Dr2rAdjoint(x, seti);
    VStarf = @(x) seti.k^2.*helmholtz2Dc2cAdjointf(x, seti);
    VqStar = @(x) conj(qROI).*VStar(x);
    VqStarf = @(x) VStarf(fftn(extendROItoCD(conj(qROI) .* restrictCDtoROI(ifftn(x),seti.ROImask),seti.ROImask)));
    
    QU = @(x) qROI .* x;
    QStarU = @(x) conj(qROI) .* x;
elseif  strcmp(seti.model, 'helmholtz3D')
    V  = @(x) seti.k^2.*helmholtz3Dr2r(x, seti); % volume potential
    Vf  = @(x) seti.k^2.*helmholtz3Dc2cf(x, seti);
    Vqf = @(x) Vf(fftn(extendROItoCD(qROI .* restrictCDtoROI(ifftn(x),seti.ROImask),seti.ROImask)));
    
    VGStar = @(x) seti.k^2.*helmholtz3Dr2dataAdjoint(x, seti);
    
    VStar  = @(x) seti.k^2.*helmholtz3Dr2rAdjoint(x, seti);
    VStarf = @(x) seti.k^2.*helmholtz3Dc2cAdjointf(x, seti);
    VqStar = @(x) conj(qROI).*VStar(x);
    VqStarf = @(x) VStarf(fftn(extendROItoCD(conj(qROI) .* restrictCDtoROI(ifftn(x),seti.ROImask),seti.ROImask)));
    
    QU = @(x) qROI .* x;
    QStarU = @(x) conj(qROI) .* x;
elseif  strcmp(seti.model, 'helmholtzHMode2D') % not supported in public version
    V  = @(x) helmholtzHMode2Dr2r(x, qROI, seti); % volume potential
    Vf  = @(x) helmholtzHMode2Dc2cf(x, qROI, seti);
    Vqf =  Vf;
    
    VGStar = @(x) helmholtzHMode2Dr2dataAdjoint(x, qROI, seti);
    
    VStar  = @(x) helmholtzHMode2Dr2rAdjoint(x, qROI, seti);
    VqStar = @(x) VStar(x);
    
    VStarf  = @(x) helmholtzHMode2Dc2cAdjointf(x, qROI, seti);
    VqStarf = @(x) VStarf(x);
    
    QU = @(x) x;
    QStarU = @(x) x;
else
    disp('intOpsFuncs.m: Error - model is not (yet ... ?) implemented');
end

%%

S = @(s) simo(qROI, s, seti); % simo has support for two-grid (two-grid not supported in public version)

if seti.mCD == 0
    TStar  = @(x) solveLippmannSchwinger(VqStar,x,seti);
else
    [~, ~, ~, ~, ~, ~, ~, ~, VqStarfM] = intOpsFuncs(seti.setiM.qROI, seti.setiM);
    
    TStar = @(x) solveLippmannSchwinger(VqStar, x, seti, 'twoGrid', VqStarf, VqStarfM);
end

end
