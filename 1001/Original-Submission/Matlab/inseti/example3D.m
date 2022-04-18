%% inseti: example 3D
% Example of input parameters to run |start| not with default parameters.
% See Also: <example.html>.

% -- dirname: suffix
seti.dirSuffix = '_test';

% -- setGrid
seti.dim = 3;
seti.rCD = 2;
seti.nCD = 256;

% -- setKernel
seti.k = 10;
seti.model = 'helmholtz';

% -- setContrast
seti.contrast = 'twoTripods3D';

% -- expSetup (set experimental set-up)
seti.incPntsType = 'sphereFibo';
seti.measPntsType = 'sphereFibo';
seti.incNb = 31;
seti.measNb = 31;
seti.radSrc = 5;
seti.radMeas = 5;
seti.incType = 'pointSource';
seti.measType = 'nearField';

% -- reconstruction
seti.invNo = 6;
seti.delta = 0.01;
seti.physBounds = [-1,3,0,3];
seti.alpha = 500;
seti.beta = 1E-5;
seti.useDis = 1;
seti.tau = 1.25;
seti.nOut = 30;

% -- reconstruction with PDA
seti.pdaN = 50;
seti.pdaStepsize = 'fix';

% -- figures and files
seti.plotFreq = 1;
seti.plotPublish = 0;

seti.usecbarlim = 0;
seti.cbarlim = [-0.2, 1.4];

seti.savepng = 1;
seti.saveepsc = 0;
seti.savefig = 0;

seti.savedata = 1;

