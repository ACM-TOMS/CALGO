%% fresnel
% Example of input parameters to run |start| using Fresnel data.
%
%% Syntax and Example
%   inseti = 'fresnel';
%   start
%
% * General information about using |inseti| are in <example.html>
% * Make sure that the experimental data from Institute Fresnel 
%   are in the correct folder, in this case: 
%   |inexpdata/fresnel_opus_1/twodielTM_8f.exp|.
% * Real-world data from Institute Fresnel are not part of this package. 
% You have to download them separately, see <readRAWData.html>.
% 
%
%% Input Arguments
% Specific input arguments to use the computational framework with data
% from Institute Fresnel:
%
% * seti.expData = 'fresnel'; : in general to use experimental data from Institute Fresnel
%                               (currently no other data is available).
% * seti.nuMax  :   2*nuMax+1 coefficients are computed for polynomial
%                   approximation, see <matchIncField.html> for details.
% * seti.ampCalc    :   method to compute coefficients, we recommend 1,
%                       see <matchIncField.html> for details.
% * seti.fresnelFreq    :   Frequency in Hz
%                           (make sure that data for this frequency is available, see [1] and [2])
% * seti.fresnelFile    :   Path to *.exp-file with data from Institute Fresnel
%
% For other fields in struct |seti| see <setiRef.html>.
%
%% References
%
% * [1] Kamal Belkebir and Marc Saillard. 
%       Special section on testing inversion algorithms against experimental data. 
%       _Inverse Problems_, 17(6):1565-1571, 2001.
% * [2] Jean-Michel Geffrin, Pierre Sabouroux, and Christelle Eyraud.
%       Free space experimental scattering database continuation: experimental set-up and measurement precision.
%       _Inverse Problems_, 21(6):S117, 2005.
%
%
%% See Also
% * <example.html>
% * <readRAWData.html>
% * <matchIncField.html>
% * <setiRef.html>
%
%% Code

% -- setGrid
seti.dim = 2;
seti.rCD = 0.2;
seti.nCD = 256;

% -- setKernel
seti.model = 'helmholtz';

% -- setExpData
seti.expData = 'fresnel';
seti.nuMax = 10;
seti.ampCalc = 1;

% -- reconstruction
seti.beta = 1E-5;
seti.alpha = 500;
seti.tau = 1.6; 

seti.physBounds = [-1,3,0,1];

seti.nOut = 20;
seti.useDis = 1;

seti.pdaStepsize = 'fix';
seti.pdaN = 10;

% -- figures and files
seti.usecbarlim = 1;
seti.cbarlim = [-0.40,2.40];

seti.plotPublish = 1;
seti.plotFreq = 10;
seti.plotFreqiPda = 0;
seti.savedata = 0;

% -- Fresnel: 2 dielectrics with 5 GHz
seti.delta = 0.25;
seti.fresnelFreq = 5*1E9;
seti.fresnelFile = 'inexpdata/fresnel_opus_1/twodielTM_8f.exp';
