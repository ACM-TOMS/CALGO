% ticGuide = tic;
init;                   % Initialization (addpath...)

seti.dim = 3;
seti.rCD = 2.0;
seti.k = 10;
seti.contrast = 'twoTripods3D';
seti.incNb = 35;
seti.measNb = 35;
seti.radSrc = 5;
seti.radMeas = 5;
seti.tau = 1.25;

seti = setData(seti);   % Set data: experimental set-up, contrast, simulate data
seti = setRecon(seti);  % Settings for the variational reconstruction
seti = recon(seti);     % Variational reconstruction (process)
% tocGuide = toc(ticGuide)

figure(1); contourPlotROI(seti.qROIexact, seti, 'real');
% set(gca,'FontSize',20); print(1,'-depsc','guideRecon3D1.eps');

figure(2); contourPlotROI(seti.qROIcomp, seti, 'real');
% set(gca,'FontSize',20); print(2,'-depsc','guideRecon3D2.eps');
% set(gca,'FontSize',20); print(2,'-dpng','guideRecon3D2dpiDefault.png');
% set(gca,'FontSize',20); print(2,'-dpng','-r600','guideRecon3D2dpi600.png');

fprintf('Reconstruction stopped after %i outer iterations...\n', seti.iOutStop);
fprintf('... with relative discrepancy %.4f...\n', seti.dis(seti.iOutStop));
fprintf('... and relative error %.4f.\n', seti.err(seti.iOutStop));

% -- Result --
% >> init; guideRecon3D
% 
% tocGuide =
% 
%    1.4730e+04 % 4.0916 h
% 
% Reconstruction stopped after 9 outer iterations...
% ... with relative discrepancy 0.0125...
% ... and relative error 0.6590.
% >> 

