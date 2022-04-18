init;
guideSetGrid;

q2D = fresnel_op1_twodielTM(seti.gridROI(1,:), seti.gridROI(2,:),seti);
x = seti.gridROI(1,:); y = seti.gridROI(2,:);
figure(1); imagesc(x,y,real(reshape(q2D,[seti.nROI seti.nROI]))); axis xy; colorbar;
% set(gca,'FontSize',20); axis square; print(1,'-depsc','guideWorkingFresnelContrast.eps');
