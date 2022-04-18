init;
setInput;
% Definition of fields of seti
seti.contrast = 'cornerBallSparseMod2D';
seti.radSrc = 5;
seti.radMeas = 6;

seti = setData(seti,4,2);
seti = setRecon(seti,4,2);
seti = recon(seti,4,2);

