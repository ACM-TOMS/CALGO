%% setData
% Sets geometry and experimental set-up, 
% simulates the scattered field and computes noisy data 
% (or loads real-world data).
%
%% Syntax
%
%   seti = setData(seti)
%   seti = setData(seti,dispDepth)
%   seti = setData(seti,dispDepth,out)
%
%% Description
%
% |seti = setData(seti)| essentially 
% sets the grids (seti.grid and seti.gridROI)
% as well as the transmitters and receivers positions (seti.incPnts and
% seti.measPnts), simulates the scattered field (seti.FmeasExact) and
% computes noisy data (seti.FmeasDelta) (or loads real-world data).
%
% |seti = setData(seti,dispDepth)| does the same as |seti = setData(seti)|
% but allows to control the depth of displayed messages by |dispDepth|.
%
% |seti = setData(seti,dispDepth,out)| does the same as |seti = setData(seti,dispDepth)| in 
% case of |out = 0|, but plots the figures in case of |out = 1| and 
% additionally saves them in case of |out = 2|.
%
% * If |seti.expData = 'fresnel'| real-world data from 
% Institute Fresnel are loaded (see <loadData.html>), i.e. the measurements of the incident 
% fields and total fields are loaded, the corresponding geometry and 
% simulation is done, and the incident field is matched.
%
% * If |seti.loadFmeas| contains a path, |FmeasExact| and |FmeasDelta| are
% loaded.
% * Note that |FmeasDelta| is only used if |seti.useFmeasDelta| was set to 1
% (otherwise data with noise is generated again).
%
% * If no field |expData| of |seti| exists and |out| equals 2, the exact
% and simulated data |seti.FmeasExact| and |seti.FmeasDelta| are saved as
% |FmeasExact| and |FmeasDelta| in a file |save_Fmeas.mat|. (Exactly this
% file can be loaded by |seti.loadFmeas|.)
%
% *Structure of the function*
%
% * The set of geometry and simulation (see <setGeomSim.html>) essentially sets
% general settings for figures, 
% the grids (computational domain (CD) and region of interest (ROI)), 
% the kernel, the experimental set-up, and the predefined contrast.
% * The the scattered field is simulated by <mimo.html>.
% * Noise is added by <addNoise.html>.
%
% 
%% Example
%
% *Example 1*
%
%  init;                    % initialization
%  seti = struct;           % create empty structural array
%  seti = setData(seti);    % compute but do not plot or store figures
%
% *Example 2*
%
%  init;
%  seti = struct;
%  seti = setData(seti,4);  % display output messages
%
% *Example 3*
%
%  init;
%  setInput;   % set input and make directories for output
%  seti = setData(seti,4,2);    % display output messages; plot and save files and figures
%
%% Input Arguments
%
% For a description of the following input arguments see <setInput.html>.
%
% * |inseti|
% * |seti|
% * |usevaralpha|
% * |usevarbeta|
% * |usevardelta|
% * |usevartol|
%
% Note that many of the output arguments can by set in the struct seti as 
% input arguments. This fields are marked by 
% "Can be set by user (otherwise default is set)".
%
% We only mention explictly the input arguments to deal with experimentally
% measured data from Institute Fresnel.
%
% * seti.expData    :   Set 'fresnel' to load real-world data
%                       (experimentally measured) from Institute Fresnel.
%
% *The following parameters can be set by user only in case of 
%  seti.expData = 'fresnel'*
% 
% For further details see <loadData.html> and <matchIncField.html>.
%
% * seti.fresnelFreq    :   frequency of Fresnel data set in Hz
%                           (default: 5e+09)
% * seti.fresnelFile    :   path to real-world data
%                           (default:
%                           'inexpdata/fresnel_opus_1/twodielTM_8f.exp')
%                           (Make sure to use data with |TM| in the filemname,
%                           which stands for transverse magnetic (TM) 
%                           polarization, because the Helmholtz equation
%                           only models electromagnetic waves in the TM
%                           case, see Section 2 in [1].)
% * seti.nuMax          :   match incident field parameter (default: 7)
% * seti.ampCalc        :   method to compute the coefficients c (1, 2 or 3),
%                           we recommend to use method 1.
%
% *Optional Input Argument*
%
% * dispDepth     : depth of displayed messages in dependence of
% |dispDepth| (number between 0 and 5).
% * out           : output depth: no figure (0), plot figure (1), 
%                   plot and save figure (2).
%
%% Output Arguments
%
% * seti    :   structure array
%
% For the following output arguments see <setGeomSim.html>:
%
% * seti.tol
%
%
% For the following output arguments see Subfunction |setFigureSettings| in
% <setGeomSim.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.plotFreq       :   default: 1
% * seti.usecbarlim     :   default: 1
% * seti.cbarlim        :   default: [-0.2 1.4]
% * seti.plotPublish    :   default: 0
% * seti.pubFontSize    :   default: 20
% * seti.plotVisible    :   default: 'off'
% * seti.savepng        :   default: 1
% * seti.saveepsc       :   default: 0
% * seti.savefig        :   default: 0
% * seti.plotFreqiPda   :   default: 0
%
% For the following output arguments see <setGrid.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.dim            :   default: 2
% * seti.rCD            :   default: 0.2
% * seti.nCD            :   default: 256
%
% _Are set automatically_:
%
% * seti.h              :   2*rCD/nCD
% * seti.dV
% * seti.grid           :   size: seti.dim x seti.nCD^seti.dim
% * seti.ballMask       :   size: seti.nCD x seti.nCD (logical)
% * seti.nROI
% * seti.gridROI        :   size: seti.dim x seti.nROI^seti.dim
% * seti.ROImask        :   size: nCD x nCD (logical)
% * seti.ballMaskROI    :   size: nROI x nROI (logical)
%
% For the following output arguments see <setReshapeVecMat.html>:
%
% _Are set automatically:_
%
% * seti.GROI
% * seti.GCD
% * seti.G
% * seti.iG
%
% For the following output arguments see <setGridScale.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.gscale     : default: 0
% * seti.nCDinv     :
%
% _Are set automatically:_
%
% * seti.nInv       :
% * seti.hInv       :
% * seti.dVinv      :
% * seti.GInv       :
% * seti.GU         :
% * seti.GD         :
%
% For the following output arguments see <setIdImagReal.html>:
% 
% _Are set automatically:_
%
% * seti.S  :   
% * seti.R  :   
% * seti.I  :   
% * seti.T  :   
%
% For the following output arguments see <setKernel.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.k      : default: 250
% * seti.model  : default: 'helmholtz2D'
%
% _Are set automatically:_
%
% * seti.kHat   : size nCD x nCD
%
% For the following output arguments see <expSetup.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.incType        :   default: 'pointSource'
% * seti.measType       :   default: 'nearField'
% * seti.incNb          :   default: 35
% * seti.measNb         :   default: 35
% * seti.radSrc         :   default: 5
% * seti.radMeas        :   default: 5
% * seti.incPntsType    :   default: 'circle'
% * seti.measPntsType   :   default: 'circle'
%
% _Are set automatically:_
%
% * seti.incPnts        :   size seti.dim x seti.incNb
% * seti.dSInc          :   
% * seti.measPnts       :   size seti.dim x seti.measNb
% * seti.dSMeas         :   
% 
% * seti.incField       :   complex matrix of size 
%                           seti.nROI^seti.dim x seti.incNb
% * seti.measKer        :   complex matrix of size 
%                           seti.measNb x seti.nROI^seti.dim
%
% For the following output arguments see <setContrast.html>:
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.contrast   :   default: 'cornerBallSparse2D'
%
% _Are set automatically:_
%
% * seti.qCDexact   :   size: seti.nCD^seti.dim x 1
% * seti.qROIexact  :   size: seti.nROI^seti.dim x 1
%
%
% For the following output argument see <setGeomSim.html>:
%
% Note that seti.mCD unequal 0 is not supported in public version.
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.mCD    :   default: 0
% 
% For the following output arguments see <addNoise.html>:
% 
% _Can be set by user (otherwise default is set)_:
%
% * seti.delta      :   default: 0.01
% * seti.whichNoise :   default: 'normal'
% * seti.seed       :   default: 0
%
% The following output arguments are set in setData (this function):
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.loadFmeas      :   path to mat-file containing variables 
%                           FmeasExact and FmeasDelta.
%                           (default: empty, i.e. '').
% * seti.useFmeasDelta  :   Set it to 1 to use FmeasDelta in loadFmeas
%                           (default: 1).
%
% _Are set automatically:_
%
% * seti.FmeasExact     :   scattered field evaluated on receivers positions 
%                           (complex matrix of size seti.measNb x seti.incNb),
%                           see also <mimo.html>.
% * seti.FmeasDelta     :   FmeasExact with noise (synthetic data), see <addNoise.html>
%                           (complex matrix of size seti.measNb x seti.incNb).
%
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <loadData.html>
% * <setGeomSim.html>
% * <mimo.html>
% * <addNoise.html>
% * <setInput.html>
%
% * <matchIncField.html>
%
% * <setGrid.html>
% * <setReshapeVecMat.html>
% * <setGridScale.html>
% * <setIdImagReal.html>
% * <setKernel.html>
% * <expSetup.html>
% * <setContrast.html>
%
%
%% Code: setData
% * geometry of experimental set-up
% * grid generation
% * simulation of forward scattering and add noise to data
% * or read and process real-world data
%

function seti = setData(seti,varargin)

if nargin == 1
    dispDepth = 0;
    out = 0; % plot no figures
elseif nargin == 2
    dispDepth = varargin{1};
    out = 0;
elseif nargin == 3
    dispDepth = varargin{1};
    out = varargin{2}; % output: out = 1 plots figures, out = 2 additionally saves them
else
    error('nargin has to be between 1 and 3.');
end

if isfield(seti,'expData') % load real-world data
    if strcmp(seti.expData,'fresnel')
        %%
        % *expData: Fresnel*
        seti = checkConsisExpData(seti,dispDepth);
        if dispDepth >= 1
            disp('   Loading fresnel data from:')
            fprintf('      %s\n',seti.fresnelFile)
            fprintf('      for frequency %d GHz',seti.fresnelFreq/1E9)
        end
        seti = loadData(seti,dispDepth,out); % setGeomSim is called in loadData
        % seti.qROIexact exact contrast in ROI as a vector (or zeroes if no exact contrast available)
        % seti.FmeasDelta contains the Fresnel data
    elseif strcmp(seti.expData,'simonetti')
        %%
        % *expData: Simonetti*
        %
        % This case is not supported in public version.
        %
        disp('loading Simonetti"s data') ;
        seti = loadSimonettiData(seti); % setGeomSim is called in loadData
        seti = setGeomSim(seti,dispDepth,out); % set geometry and simulation
        
        seti.FmeasDeltaSim = mimo(seti, seti.qROIexact, 'simo');
        
        seti.FmeasDelta = seti.FmeasDelta/norm(seti.FmeasDelta)*norm(seti.FmeasDeltaSim);
        seti.FmeasExact = seti.FmeasExact/norm(seti.FmeasExact)*norm(seti.FmeasDeltaSim);

        if ~isfield(seti,'gradMat')
            N = seti.nCD;
            N2 = ceil((seti.nCD-1)/2);
            seti.gradMat = circshift(diag((0:N-1)-N2), -N2*[1, 1]);
        end
    else
        error('seti.expData not found .. ?!?');
    end
else
    %%
    % *Set data structures for geometry and simulation*
    %
    if dispDepth >= 1
        disp('-- setGeomSim --')
        disp(' ')
    end
    seti = setGeomSim(seti,dispDepth,out); % set geometry and simulation

    % Create gradient structures if seti.model='helmholtzHMode2D'
    % Note that ...HMode... is not supported in public version.
    if strcmp(seti.model,'helmholtzHMode2D')
        if ~isfield(seti,'gradMat')
            N = seti.nCD;
            N2 = ceil((seti.nCD-1)/2);
            seti.gradMat = circshift(diag((0:N-1)-N2), -N2*[1, 1]);
        end
    end
    
    seti = checkfield(seti,'loadFmeas','',dispDepth); % default: '', so no load...
   
    if ~isempty(seti.loadFmeas)
% -- load old: start --
% 		load(sprintf('%s',seti.loadFmeas)); % loads: FmeasExact and FmeasDelta
% 		seti.FmeasExact = FmeasExact; % exact data (just loaded)
% 		clear FmeasExact k directions;
% -- load old: end --
% -- load new start --
		sloaded = load(sprintf('%s',seti.loadFmeas)); % loads: FmeasExact and FmeasDelta (in struct sloaded)
		seti.FmeasExact = sloaded.FmeasExact; % exact data (just loaded)
% -- load new end --
        if dispDepth >= 1
            disp(' - FmeasExact and FmeasDelta loaded')
        end
    else
        if dispDepth >= 1
            disp(' - FmeasExact (exact data) computation')
        end
        % Note that mimo is expensive if qROI does not have many entries with 0.
        % Define qROI before operator FF and then compute FmeasExact = FF(q).
        % Be careful with mimo and wavelets because mimo is FF(q)-FmeasDelta(!).
        ticExact = tic;
		seti.FmeasExact = mimo(seti, seti.qROIexact, 'simo');
        tocExact = toc(ticExact);
        if dispDepth >= 1
            fprintf('   Elapsed time of computation of exact data is %05.1f min.\n',tocExact/60);
        end
    end
    
    seti = checkfield(seti,'useFmeasDelta',1,dispDepth);
    
    if ~isempty(seti.loadFmeas) && seti.useFmeasDelta == 1
        if dispDepth >= 1
            disp('- use FmeasDelta (use loaded data with noise)')
        end
        seti.FmeasDelta = sloaded.FmeasDelta; % noisy data was loaded
    else
        % -- add noise to data
        if dispDepth >= 1
            disp(' - add noise to data')
        end
        [seti, seti.FmeasDelta] = addNoise(seti,seti.FmeasExact,dispDepth);
    end
    % FmeasExact and FmeasDelta are saved in start.m (first construct seti.dirname)
end

%% 
% *Save FmeasExact and FmeasDelta*
%
if ~isfield(seti,'expData') && out == 2
    saveDes = sprintf('%s/save_Fmeas.mat',seti.dirname);
    FmeasExact = seti.FmeasExact; %#ok (suppress MATLAB warning, because exact data will be saved)
    FmeasDelta = seti.FmeasDelta; %#ok (suppress MATLAB warning, because noisy data will be saved)
    save(saveDes,'FmeasExact','FmeasDelta');
    clear saveDes FmeasExact FmeasDelta;
end

end
