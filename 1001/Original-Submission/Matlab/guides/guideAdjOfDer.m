init;
seti.contrast = 'corner2D';
seti.rotation = 20;
seti = setGeomSim(seti);
qROI = seti.qROIexact;
FmeasDelta = zeros(seti.measNb,seti.incNb);
[ADFFq,seti] = adjOfDer(seti,qROI,FmeasDelta);
figure(1); imagesc(real(seti.G(ADFFq))); axis xy; colorbar;
% set(gca,'FontSize',20); axis square; print(1,'-depsc','guideAdjOfDer.eps');
