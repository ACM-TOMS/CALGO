%% Example for transmitters on a sphere (method: Fibonacci lattice)

init;
seti.dim = 3;
seti.incType = 'pointSource';
seti.incPntsType = 'sphereFibo';
seti.radSrc = 5;
seti.incNb = 50;
seti = setIncPnts(seti);
figure(1); [sx,sy,sz] = sphere; r = seti.radSrc; surf(sx*r,sy*r,sz*r); hold on;
Pnts = seti.incPnts; s = 100; scatter3(Pnts(1,:),Pnts(2,:),Pnts(3,:),s,'filled','b'); hold off; % s: size
axis square; set(gca,'FontSize',20); print(1,'-depsc','guideExpSetupTrans3D.eps');

