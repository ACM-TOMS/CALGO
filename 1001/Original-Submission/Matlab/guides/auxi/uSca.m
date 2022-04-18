function [uTotRX, uIncRX, uScaRX] = uSca(uTotRX,uIncRX,frequencyId)
% Compute the scattered field uScaRX at receivers positions:

% Conjugate the fields to adapt from time dependence exp(iwt) to exp(-iwt):
uTotRX   = conj(uTotRX);
uIncRX   = conj(uIncRX);

% FmeasDelta = uScaRX = uTotRX - uIncRX:
uTotRX = uTotRX(:,:,frequencyId);
uIncRX = uIncRX(:,:,frequencyId);

%
uScaRX = uTotRX - uIncRX; % uScaRX = uTotRX - uIncRX
end