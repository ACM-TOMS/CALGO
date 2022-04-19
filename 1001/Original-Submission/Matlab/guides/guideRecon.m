init;                   % Initialization (addpath...)
seti = struct;
seti = setData(seti);   % Set data: experimental set-up, contrast, simulate data
seti = setRecon(seti);  % Settings for the variational reconstruction
seti = recon(seti);     % Variational reconstruction (process)

figure(1); imagesc(real(seti.G(seti.qROIexact))); axis xy; colorbar;
%set(gca,'FontSize',20); axis square; caxis([-0.15 +1.42]); print(1,'-depsc','guideRecon1.eps');

figure(2); imagesc(real(seti.G(seti.qROIcomp))); axis xy; colorbar;
%set(gca,'FontSize',20); axis square; caxis([-0.15 +1.42]); print(2,'-depsc','guideRecon2.eps');

fprintf('Reconstruction stopped after %i outer iterations...\n', seti.iOutStop);
fprintf('... with relative discrepancy %.4f...\n', seti.dis(seti.iOutStop));
fprintf('... and relative error %.4f.\n', seti.err(seti.iOutStop));
