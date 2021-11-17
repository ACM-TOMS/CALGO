% Experimental Setup
% clear all;

init;

guideExpSetupTrans; % set transmitters positions
guideExpSetupRece;  % set receivers positions
guideSetGrid; % set a grid

seti.model = 'helmholtz2D';
seti.k = 250; % wave number

seti = setIncField(seti); % Evaluate incident fields on ROI
seti = setMeasKer(seti); % set up measurement kernels (such that measurement = k^2 * kernel * solution * voxelVolume)

figure(1); imagesc(real(reshape(seti.incField(:,1),[seti.nROI seti.nROI]))); axis xy; colorbar; % incident field in ROI of 1st transmitter
% set(gca,'FontSize',20); axis square; print(1,'-depsc','guideExpSetup1.eps');

figure(2); imagesc(real(reshape(seti.measKer(5,:),[seti.nROI seti.nROI]))); axis xy; colorbar; % measurement kernel on ROI of 5th receiver
% set(gca,'FontSize',20); axis square; print(2,'-depsc','guideExpSetup2.eps');

