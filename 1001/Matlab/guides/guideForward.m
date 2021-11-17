init;
seti.contrast = 'corner2D';
seti.rotation = 20;
seti = setGeomSim(seti);
figure(1); imagesc(real(seti.G(seti.qROIexact))); axis xy; colorbar;
% set(gca,'FontSize',20); axis square; print(1,'-depsc','guideForward1.eps');
[FFqMeas,FFqROI,seti] = forward(seti,seti.qROIexact);
figure(2); imagesc(real(seti.G(FFqROI(:,7)))); axis xy; colorbar; % 7-th transmitter
% set(gca,'FontSize',20); axis square; print(2,'-depsc','guideForward2.eps');
figure(3); imagesc(real(FFqMeas)); axis xy; colorbar;
% set(gca,'FontSize',20); axis square; print(3,'-depsc','guideForward3.eps');

%% add noise...
seti.delta = 0.05;
seti.whichNoise = 'normal';
seti.seed = 10;
[seti, FmeasDelta] = addNoise(seti, FFqMeas);
figure(4); imagesc(real(FFqMeas)); axis xy; colorbar;
set(gca,'FontSize',20); axis square; print(4,'-depsc','guideForward4.eps');

%% Alternative
% init;
% % Set a grid in 2D
% seti.qBall = 0.8;
% seti.rBall = 0.05;
% guideSetGrid;
% % Contrast
% qROI = referenceBall2D(seti.gridROI(1,:), seti.gridROI(2,:), seti);
% % Forward
% [FFqMeas,FFqROI,seti] = forward(seti,reshape(qROI,[seti.nROI^seti.dim 1]));
