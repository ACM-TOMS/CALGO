init;
seti.contrast = 'corner2D';
seti.rotation = 20;
seti = setGeomSim(seti);
qROI = seti.qROIexact;
[JA,JB] = derivative(seti,qROI);

DFFq = @(h) JA*diag(h)*JB;
h = ones(size(qROI))+1i*ones(size(qROI));
DFFqh = DFFq(h); % Evaluate F'(q)[h].
figure(1); imagesc(real(DFFqh)); axis xy; colorbar;
% set(gca,'FontSize',20); axis square; print(1,'-depsc','guideDerivative.eps');
