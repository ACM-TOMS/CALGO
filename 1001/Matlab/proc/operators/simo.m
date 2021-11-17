%% simo
% simo: single input and multiple output.
%
%% Syntax
%
%   [uScattRX, uScattROI] = simo(qROI, uIncROI, seti)
%
%% Description
%
% |[uScattRX, uScattROI] = simo(qROI, uIncROI, seti)| computes 
% one scattered field |uScattRX| at receivers and |uScattROI| in ROI 
% for the contrast |qROI| and one incident field |uIncROI|.
% 
%% Input Arguments
%
% * qROI    : contrast
% * uIncROI : incident field
% * seti    : structural array
%
% See also <mimo.html> for further information.
%
%% Output Arguments
%
% * uScattRX  : scattered field at RX (receivers points)
% * uScattROI : scattered field in ROI
%
% See also <mimo.html> for further information.
%
%% More About
%
% * Function maps column vectors qROI, uIncROI to column vectors
% uScattRX, uScattROI.
% * Note: simo can deal with two-grid (two-grid not available in public version).
% * Computation is described in <solveLippmannSchwinger.html>.
%
%% See Also
%
% * <mimo.html>
% * <intOpsFuncs.html>
% * <solveLippmannSchwinger.html>
%
%% Code

function [uScattRX, uScattROI] = simo(qROI, uIncROI, seti)

%%
% *Compute the scattered field on ROI: |uScattROI|*

[V, ~, ~, ~, ~, QU, ~, Vqf, ~] = intOpsFuncs(qROI, seti);

if seti.mCD == 0 % in public version seti.mCD = 0 because twoGrid is not supported
    uScattROI = solveLippmannSchwinger(@(x) V(QU(x)), V(QU(uIncROI)), seti); 
else % twoGrid
    [~, ~, ~, ~, ~, ~, ~, VqfM, ~] = intOpsFuncs(seti.setiM.qROI, seti.setiM);
    uScattROI = solveLippmannSchwinger(@(x) V(QU(x)), V(QU(uIncROI)), seti, 'twoGrid', Vqf, VqfM);
end

%%
% *Compute the scattered field at receivers positions RX: |uScattRX|*

fROI = QU(uIncROI+uScattROI);
if strcmp(seti.model, 'helmholtz2D')
    uScattRX = seti.k^2.*helmholtz2Dr2data(fROI, seti);
elseif strcmp(seti.model, 'helmholtz3D')
    uScattRX = seti.k^2.*helmholtz3Dr2data(fROI, seti);
elseif strcmp(seti.model, 'helmholtzHMode2D') % not supported in public version
    uScattRX = helmholtzHMode2Dr2data(fROI, qROI, seti);
else
    disp('simo.m: Error - model is not (yet ... ?) implemented');
end

end
