% Arrange 4+1 points on each edge of a square with length 3. Total number of points is 4*4 = 16.

init;
seti.dim = 2;
seti.measType = 'nearField';
seti.measPntsType = 'square';
seti.measNbEdge = 4;
seti.measEdgeLength = 0.8;
seti = setMeasPnts(seti);
figure(1); plot(seti.measPnts(1,:),seti.measPnts(2,:),'rs','MarkerSize',7,'MarkerFaceColor','red');

% seti.measNb
% seti.measPnts

% set(gca,'FontSize',20); axis square; xlim([-0.5 0.5]); ylim([-0.5 0.5]); print(1,'-depsc','guideExpSetupRece.eps');
