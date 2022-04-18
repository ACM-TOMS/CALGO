%% readRAWData
% Reads the raw data from Institute Fresnel (1st and 2nd opus).
%
%% Syntax
%
%   [uTotRX, uIncRX, frequencies, rTX, nTX, rRX, nRX] = readRAWData(filename)
%
%% Description
% |[uTotRX, uIncRX, frequencies, rTX, nTX, rRX, nRX] = readRAWData(filename)|
%
% * Reads the raw data from |filename|, which is the path to one of the
% experimental data files from Institute Fresnel (1st or 2nd opus).
% * Serves the total fields |uTotRX| and incident fields |uIncRX| at
% receivers positions for each transmitter for each frequency.
% * Provides the |frequencies|, the number of transmitters and receivers 
%   (|nTX| and |nRX|) as well as the radius of circles they are arranged on
%   (|rTX| and |rRX|).
%
%
%% How to get experimental data from Institute Fresnel
%
% Overview of available experimental data is given in:
%
% <http://www.fresnel.fr/3Ddatabase/database.php> (Accessed: 20160921).
% 
% *First opus*
%
% The corresponding article is [1].
%
% # Open <http://dx.doi.org/10.1088/0266-5611/17/6/301> (Accessed: 20160921).
% # Klick on "Supplementary Data".
% # Download the experimental data: |*.exp|-files.
% # Save them in |inexpdata/fresnel_opus_1|.
%
% *Second opus*
%
% The corresponding article is [2].
%
% # Open <http://dx.doi.org/10.1088/0266-5611/21/6/S09> (Accessed: 20160921).
% # Klick on "Supplementary Data".
% # Download the experimental data: |*.exp|-files.
% # Save them in |inexpdata/fresnel_opus_2|.
% 
%
%% Example
%
% *Make sure that the file |twodielTM_8f.exp| is in the folder
% |inexpdata/fresnel_opus_1|*. 
% 
% A description how to do this is given in the
% section "How to get..." in the part "First opus" above.
%
% Make sure that you are currently in the folder |proc/expData|.
%
%   filename = '../../inexpdata/fresnel_opus_1/twodielTM_8f.exp';
%   [uTotRX, uIncRX, frequencies, rTX, nTX, rRX, nRX] = readRAWData(filename);
%
% This example reads the experimental data from |twodielTM_8f.exp|.
%
% It provides experimental data of two dielectric cylinders.
% See <incontrastsRef.html> and [1] for details.
%
%
%% Input Arguments
%
% * filename    :   path to file with experimental data from Institute Fresnel
%
%
%% Output Arguments
%
% * |uTotRX|    :   *total field* (i.e. with obstacle) at receivers positions 
%                   for each transmitter for each frequency
%                   (complex array of size nRX x nTX x _number of frequencies_).
% * |uIncRX|    :   *incident field* (i.e. without obstacle) at receivers positions 
%                   for each transmitter for each frequency
%                   (complex array of size nRX x nTX x _number of frequencies_).
% * |frequencies|   :   Available frequencies as a vector of size N x 1.
% * |nTX|           :   number of transmitters
% * |nRX|           :   number of receivers 
% * |rTX|           :   radius of circle transmitters are arranged on in meters
% * |rRX|           :   radius of circle receivers are arranged on in meters
%
%
%% More About
%
% * For *metrologically reasons* it is not possible to measure at some
% receiver positions depending on the active transmitter.
% * This *entries* in |uTotRX| and |uIncRX| are marked with *NaN* (not a number).
% * If you use the data in the provided computational framework, make sure
% to use only transverse magnetic (TM) polarization and not TE
% polarization, because for electromagnetic waves the Helmholtz equation
% does only model the TM case, see Section 2 in [3].
%
% *Fixed settings for Fresnel data*
%
% * seti.dim = 2;                   : opus 1 and 2 are two dimensional
% * seti.incType = 'pointSource';   : transmitters are point sources
% * seti.measType = 'nearField';    : receivers measure near field data
%
% See also <fresnel.html>, <checkConsisExpData.html>
%
%% References
%
% * [1] Kamal Belkebir and Marc Saillard. 
%       Special section on testing inversion algorithms against experimental data. 
%       _Inverse Problems_, 17(6):1565-1571, 2001.
% * [2] Jean-Michel Geffrin, Pierre Sabouroux, and Christelle Eyraud.
%       Free space experimental scattering database continuation: experimental set-up and measurement precision.
%       _Inverse Problems_, 21(6):S117, 2005.
% * [3] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%
%% See Also
% * <incontrastsRef.html>
% * <fresnel.html>
% * <loadData.html>
% * <checkConsisExpData.html>
%
%% Code
function [uTotRX, uIncRX, frequencies, rTX, nTX, rRX, nRX] = readRAWData(filename)

dataType = whichData(filename);
rawData = importdata(filename, ' ', 10);
rawData = rawData.data;

% Distances in meters
if regexpi(dataType, 'fresnel_opus_1')
    rTX = 0.72;
    rRX = 0.76;
end
if regexpi(dataType, 'fresnel_opus_2_45deg|fresnel_opus_2_20deg')
    rTX = 1.67;
    rRX = 1.67;
end

% Number of transmitters and receivers  (angular steps in degrees)
if regexpi(dataType, 'fresnel_opus_1') 
    nTX = 36; % (10 degrees)
    nRX = 72; % (5 degrees)
end
if regexpi(dataType, 'fresnel_opus_2_45deg')
    nTX = 8; % (45 degrees)
    nRX = 360; % (1 degree)
end
if regexpi(dataType, 'fresnel_opus_2_20deg')
    nTX = 18; % (20 degrees)
    nRX = 360; % (1 degree)
end

% Frequencies (in Hz) (convert GHz in Hz)
frequencies = sort(unique(rawData(:,3)))*1E9;

% Total and incident field at reveivers points
nf = length(frequencies);
uTotRX = NaN(nRX,nTX,nf); % total field
uIncRX = NaN(nRX,nTX,nf); % incident field

for i = 1:size(rawData,1)
    t = rawData(i,:);
    iTX   = t(1);
    iRX   = t(2);
    iFreq = find(frequencies == t(3)*1E9);

    uTotRX(iRX,iTX,iFreq) =  t(4)+1i*t(5);
    uIncRX(iRX,iTX,iFreq) =  t(6)+1i*t(7);
end

end

function type = whichData(filename)
if ~isempty(regexpi(filename,strcat('dielTM_dec4f.exp|dielTM_dec8f.exp|',...
        'rectTE_8f.exp|rectTM_cent.exp|rectTM_dece.exp|twodielTM_4f.exp|',...
        'twodielTM_8f.exp|uTM_shaped.exp')))
    type = 'fresnel_opus_1';
end
if ~isempty(regexpi(filename,strcat('FoamDielExtTE.exp|FoamDielExtTM.exp|',...
        'FoamDielIntTE.exp|FoamDielIntTM.exp|')))
    type = 'fresnel_opus_2_45deg';
end
if ~isempty(regexpi(filename,strcat('FoamMetExtTE.exp|',...
        'FoamMetExtTM.exp|FoamTwinDielTE.exp|FoamTwinDielTM.exp')))
    type = 'fresnel_opus_2_20deg';
end
end
