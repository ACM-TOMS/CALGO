init;
seti.qBall = 0.8;
seti.rBall = 0.05;
seti.contrast = 'referenceBall2D';
seti = setGeomSim(seti);

imagesc(real(seti.G(seti.qROIexact))); axis xy; colorbar;
