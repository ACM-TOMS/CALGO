%% loadData
% Loads real-world (experimentally measured) data from Institute Fresnel. 
% (A consistent name would be "expData" instead of "loadData").
%
%% Syntax
%   seti = loadData(seti)
%   seti = loadData(seti,dispDepth)
%   seti = loadData(seti,dispDepth,out)
%
%% Description
% |seti = loadData(seti)| loads data from the Fresnel database via 
% function |readRAWData| and stores them in structure array |seti|.
%
% |seti = loadData(seti,dispDepth)| does the same but allows to control the
% displayed messages by |dispDepth|.
%
% |seti = loadData(seti,dispDepth,out)| does the same but transfers the argument
% |out| for output depth (plot no figure (0), plot figures (1), plot and
% save figure (2)).
%
%% Example
%
% Make sure that the database from 1st opus of Institute Fresnel is available
% as described in <readRAWData.html>.
%
%   init;
%   seti.expData = 'fresnel';
%   seti = checkConsisExpData(seti);
%   seti = loadData(seti);
%
%% Input Arguments
%
% * |seti|      : structure array
%
% The input arguments are the input and output arguments of
% |checkConsisExpData| in case of |seti.expData = 'fresnel'|, i.e.:
%
% * |seti.rCD|  : Size of computational domain [-rCD,rCD)^dim
%                 (default: 0.2 m), see <setGrid.html>.
% * |seti.nCD|  : Number of discretization points for each dimension of CD
%                 (default: 256), see <setGrid.html>.
% * |seti.fresnelFreq|  : Frequency of the Fresnel data set
%                         (default: 5 GHz = |5*1E9|)
% * |seti.fresnelFile|  : path to the file with data from Institute Fresnel
%                         (default: 'inexpdata/fresnel_opus_1/twodielTM_8f.exp').
%
% * |seti.nuMax|    : match incident field using Hankel functions of first kind and orders
%                     $\nu = -\texttt{nuMax}, ..., -1, 0, 1, ...,
%                     \texttt{nuMax}$ (default: 7)
% * |seti.ampCalc|  : method to compute the coefficients c (1, 2 or 3),
%                     (default: 1), see <matchIncField.html>.
%
% * |seti.dim = 2|  : The dimension of the problem is 2.
% * |seti.incType =  'pointSource'| : The type of incident field is set to point sources.
% * |seti.measType = 'nearField'|   : The measurement type is set to near field data.
%
% *Optional Input Arguments*
%
% * dispDepth   : depth of displayed messages (number between 0 and 5).
% * out         : depth of output: no figure (0), plot figure (1), plot and save figure (2).
%
%% Output Arguments
% Most important output arguments, i.e. fields in structure array |seti|.
%
% * |seti.incNb|      :   number of transmitters (number of incident fields)
% * |seti.measNb|     :   number of receivers (number of measurements)
% * |seti.radSrc|     :   radius of circle transmitters are arranged on
% * |seti.radMeas|    :   radius of circle receivers are arranged on
%
% * |seti.k|          :   wave number ($k = 2 \pi f/c$ with frequency f and light velocity c)
% * |seti.FmeasDelta| :   scattered field at receivers' positions (i.e. 'the data') for each transmitter
%                         (uScaRX = uTotRX - uIncRX, i.e. total field minus incident field)
%                         (complex matrix of size seti.measNb x seti.incNb)
% * |seti.incField|   :   Incident field in the region of interest (ROI)
%                         for each transmitter
%                         (complex matrix of size seti.nROI^seti.dim x seti.incNb)
%
%% More About
% This function does the following steps:
%
% # Read experimental data from Institute Fresnel.
% # Choose data to specific frequency.
% # Compute wave number.
% # Compute the scattered field |FmeasDelta| at receivers positions.
% # Set geometry and simulation to loaded data.
% # Compute the incident field |incField| on region of interest (ROI).
%
% Some more comments are in the code.
%
%% See Also
% * <expData.html>
% * <readRAWData.html>
% * <setGeomSim.html>
% * <matchIncField.html>

%% Code: loadData
function seti = loadData(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
    out = 0;
elseif nargin == 3
    dispDepth = varargin{1};
    out = varargin{2};
else
    dispDepth = 0;
    out = 0;
end

%%
% *Read experimental data from Institute Fresnel*
%
% * uTotRX and uIncRX contain NaN-values, when there are no Fresnel data.
% * Used transmitter and receiver position IDs will be marked
%   with activeTX and activeRX (note that all TX are active, but not all RX at once)

filename = seti.fresnelFile;    % filename of the data
frequencyHz = seti.fresnelFreq; % frequency in Hz
[uTotRX, uIncRX, frequencies, seti.radSrc, seti.incNb, seti.radMeas, seti.measNb] = readRAWData(filename);
% frequencies : contains all possible frequencies in Hz

%%
% *Choose data to specific frequency*
frequencyId = find(frequencies == frequencyHz);
if dispDepth >= 1
    fprintf('   FrequencyHz = %g\n',frequencyHz);
    fprintf('   frequencyId = %g\n',frequencyId);
end

if numel(frequencyId) == 0
    frequencyId = numel(frequencies);
    frequencyHz = frequencies(frequencyId);
    if dispDepth >= 1
        fprintf('   Desired frequency not found, using %1.0f GHz instead.\n', frequencyHz/1E9);
    end
end

%%
% *Compute wave number*

%c = physconst('LightSpeed'); % To use 'physconst', you might need: physconst - Phased Array System Toolbox
cLight = 2.99792458E8; % light velocity in vacuum
seti.k = 2*pi*frequencyHz/cLight;

%%
% *Compute the scattered field |FmeasDelta| at receivers positions*
%
% * Conjugate the fields to adapt from time dependence exp(iwt) to exp(-iwt).
% * FmeasDelta = uScaRX = uTotRX - uIncRX

% Conjugate the fields to adapt from time dependence exp(iwt) to exp(-iwt).
uTotRX   = conj(uTotRX);
uIncRX   = conj(uIncRX);

uTotRX = uTotRX(:,:,frequencyId);
uIncRX = uIncRX(:,:,frequencyId);

seti.FmeasDelta = uTotRX - uIncRX; % uScaRX = uTotRX - uIncRX

% Workaround to deal with undefined measurements:
% undefined values are set to zero.
seti.FmeasDelta(isnan(seti.FmeasDelta)) = 0;

%%
% *Set geometry and simulation to loaded data*
%
% See also <setGeomSim.html>

disp(' - Set geometry and simulation to loaded data')
seti = setGeomSim(seti,dispDepth,out);

%%
% *Compute the incident field |incField| on region of interest (ROI)*
%
% See also <matchIncField.html>.
%
[uIncROI,errC] = matchIncField(uIncRX,seti,'ROI');
seti.errC = errC;
seti.incField = 1/seti.dSInc*uIncROI;
% This works only in case of equidistant transmitters because 1/seti.dSInc.
% The factor 1/seti.dSInc is important because in mimo.m:
% SL = dSInc(jj)*seti.incField(:,jj);


%% Code: Compare methods to compute incField on ROI (some unused tests...)

% incField on CD is to compare
if 0

    if 0
        [uIncROI,~] = matchIncField(uIncRX,seti,'ROI');
        incFieldROI = uIncROI;
        incFieldROIAltern = matchIncidentFieldAlternative(uIncRX,seti,'ROI');
        figure(101); imagesc(real(seti.G(incFieldROI(:,1)))); title('incField on ROI, inc 1'); colorbar; axis xy;
        figure(102); imagesc(real(seti.G(incFieldROIAltern(:,1)))); title('incField on ROI, inc 1, alternative'); colorbar; axis xy;
    end
    
    % different: a small shift and other scale factor...
    % use grid and nCD to have a better comparison with Sources...

    % using nCD (only to plot it, not used later)
    [uIncCD,errC] = matchIncField(uIncRX,seti,'CD');
    seti.incFieldCD = uIncCD;
    %figure(103); imagesc(real(seti.GCD(incFieldCD(:,1)))); title('real part of incField on CD, inc 1'); colorbar; axis xy;
    if 0
        incFieldCDAltern = matchIncidentFieldAlternative(uIncRX,seti,'CD');
        figure(104); imagesc(real(seti.GCD(incFieldCDAltern(:,1)))); axis xy; title('real part of incField on CD, inc 1, alternative'); colorbar
        figure(105); imagesc(imag(seti.GCD(incFieldCD(:,1)))); axis xy; title('imag part of incField on CD, inc 1'); colorbar
        figure(106); imagesc(imag(seti.GCD(incFieldCDAltern(:,1)))); axis xy; title('imag part of incField on CD, inc 1, alternative'); colorbar
        ifield = incFieldCD(:,1);
        [phi,~] = cart2pol(real(ifield), imag(ifield));
        ifield = incFieldCDAltern(:,1);
        [phiAltern,~] = cart2pol(real(ifield), imag(ifield));
        figure(107); imagesc(seti.GCD(phi)); axis xy; title('phase of incField on CD, inc 1'); colorbar
        figure(108); imagesc(seti.GCD(phiAltern)); axis xy; title('phase of incField on CD, inc 1, alternative'); colorbar
    end

end

end
