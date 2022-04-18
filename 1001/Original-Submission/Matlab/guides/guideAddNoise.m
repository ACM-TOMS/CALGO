init;
FmeasExact = rand(10,5) + 1i*rand(10,5);
seti.dSMeas = 1;

seti.delta = 0.05;
seti.whichNoise = 'normal';
seti.seed = 10;

[seti, FmeasDelta] = addNoise(seti, FmeasExact);
