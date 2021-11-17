function [seti,uTotRX,uIncRX,uScaRX] = matchIncFieldTrans(filename,frequencyHz)

% Read Fresnel data
%filename = 'inexpdata/fresnel_opus_1/twodielTM_8f.exp';
[uTotRX, uIncRX, frequencies, rTX, nTX, rRX, nRX] = readRAWData(filename);

[frequencyId,frequencyHz] = freqID(frequencies,frequencyHz); % Choose data to specific frequency, e.g. 5 GHz via 5*1E9
seti.k = freqToWaveNumber(frequencyHz); % Compute wave number

% Choose fields at specific frequency and conjugate the fields to adapt
% from time dependence exp(iwt) to exp(-iwt):
[uTotRX, uIncRX, uScaRX] = uSca(uTotRX,uIncRX,frequencyId);

% Transfer experimental set-up into structural array seti
seti.radSrc = rTX;
seti.incNb = nTX;
seti.radMeas = rRX;
seti.measNb = nRX;

% Experimental Setup
seti.dim = 2;               % two dimensional problem
% Set transmitters positions:
seti.incType = 'pointSource';
seti.incPntsType = 'circle';
seti = setIncPnts(seti);    % seti.incPnts contains the coordinates
% Set receivers positions:
seti.measType = 'nearField';
seti.measPntsType = 'circle';
seti = setMeasPnts(seti);   % seti.measPnts contains the coordinates

% Set grid
seti.rCD = 0.2;         % size of computational domain [-rCD,rCD)^dim
seti.nCD = 256;         % number of discretization points for each dimension of CD
seti = setGrid(seti);   % grid is stored in seti.gridROI and discretization points in each dimension in seti.nROI

end