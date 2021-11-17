init;

seti.rBall = 0.05;
seti.qBall = 0.8;

guideSetGrid;
% Evaluate predefined contrast referenceBall2D in the grid
q2D = referenceBall2D(seti.gridROI(1,:), seti.gridROI(2,:), seti);
x = seti.gridROI(1,:); y = seti.gridROI(2,:);
figure(1); imagesc(x,y,real(reshape(q2D,[seti.nROI seti.nROI]))); axis xy; colorbar;
% set(gca,'FontSize',20); axis square; print(1,'-depsc','guideContrast2D.eps');

guideSetGrid3D;
q3D = referenceBall3D(seti.gridROI(1,:), seti.gridROI(2,:), seti.gridROI(3,:), seti);
figure(2); contourPlotROI(q3D, seti, 'real');
% set(gca,'FontSize',20); axis square; print(2,'-depsc','guideContrast3D.eps');
