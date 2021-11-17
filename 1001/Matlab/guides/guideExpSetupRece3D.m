
%% Example for receivers on a sphere (method: latitude-longitude lattice)

init;
seti.dim = 3;
seti.measType = 'nearField';
seti.measPntsType = 'sphereLatLon';
seti.radMeas = 5;
seti.measNb = 50;
seti = setMeasPnts(seti);
figure(1); [sx,sy,sz] = sphere; r = seti.radMeas; surf(sx*r,sy*r,sz*r); hold on;
Pnts = seti.measPnts; s = 100; scatter3(Pnts(1,:),Pnts(2,:),Pnts(3,:),s,'filled','rs'); hold off;
axis square; set(gca,'FontSize',20); print(1,'-depsc','guideExpSetupRece3D.eps');

