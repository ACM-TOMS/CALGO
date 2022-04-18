init;                   % Initialization (addpath...)
seti.useTolOut = 1;     % Use outer tolerance principle to stop the inner iteration
seti = setData(seti);   % Set data: experimental set-up, contrast, simulate data
seti = setRecon(seti);  % Settings for the variational reconstruction
seti = recon(seti);     % Variational reconstruction (process)

imagesc(real(seti.G(seti.qROIcomp))); axis xy; colorbar;
