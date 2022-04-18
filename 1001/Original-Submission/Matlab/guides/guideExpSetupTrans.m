% Set transmitters positions:
init;
seti.dim = 2;
seti.incType = 'pointSource';
seti.incPntsType = 'circle';
seti.radSrc = 0.2;
seti.incNb = 12;
seti = setIncPnts(seti);    % seti.incPnts contains the coordinates
figure(1); plot(seti.incPnts(1,:),seti.incPnts(2,:),'b.','MarkerSize',30);

% set(gca,'FontSize',20); axis square; xlim([-0.3 0.3]); ylim([-0.3 0.3]); print(1,'-depsc','guideExpSetupTrans.eps');
