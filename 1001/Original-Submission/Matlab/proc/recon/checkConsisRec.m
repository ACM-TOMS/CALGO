%% checkConsisRec
%
% The input parameters for the variational reconstruction process 
% are checked for consistency and set or corrected whenever necessary.
%
%% Syntax
%
%   seti = checkConsisRec(seti)
%   seti = checkConsisRec(seti,dispDepth)
%
%% Description
%
% |seti = checkConsisRec(seti)| checks the existence and consistency of
% some field in structure array seti.
%
% |seti = checkConsisRec(seti,dispDepth)| does the same but allows to
% control the depth of displayed messages by |dispDepth| (0: no, 1 or greater: yes).
%
% See the section Code to get links to parameters' references.
%
%% Input Arguments
%
% * |seti|        :   struct seti
%
% *Optional Input Argument*
%
% * |dispDepth|   : Depth of displayed messages (0: no, 1 or greater: yes).
%
% *Optional Input Arguments of structure |seti|*
%
% If the fields does not exist, default values are set in this function.
%
% * |seti.useWavelet|   : Use wavelets for reconstruction (default: 0).
%                         Note that wavelets are not supported in the 
%                         public version of this code. Therefore the
%                         wavelet-specific fields |wavWeight|,
%                         |genWhiteNoiseNew|, |adaptWeightsMan|, |wavelet|,
%                         |extension|, |smin| and |samples| are not
%                         documented.
% * |seti.recname|      : What is reconstructed? A |'constrast'| (default) 
%                         or a |'source'|. Influences the title in
%                         <subplots.html>.
% * |seti.physBounds|   : Bounds for real/imaginary part of contrast: 
%                         [reMin, reMax, imMin, imMax] 
%                         (default: [-1,3,0,3])
%                         (consider also the penalty for physical bounds in
%                         <pda.html> Section "More About".)
% * |seti.tolDelta|     : Tolerance for physical bounds' delta-function
%                         (default: 1E-14) in <setFuncsPda.html>.
% * |seti.useDis|       : If it is 1 (default), the discrepancy principle
%                         is used to stop the outer iteration. (Otherwise 0).
% * |seti.nOut|         : Maximal number of outer iterations (default: 30).
%
% If the primal-dual algorithm is used (|seti.inv = 'pda')
% and the following fields does not exist, default values are set in this function.
%
% * |seti.pdaN|         : Number of inner iterations (PDA) (default: 50) (see <pda.html>).
% * |seti.pdaStepsize|  : Method to choose primal and dual stepsizes
%                         (|'fix'| (default) or |'adaptive'|)
%                         see <pdaChoosingStepsizes.html>.
% * |seti.vartheta|     : Parameter in (0,1) used in case of adaptive stepsizes,
%                         see <pdaChoosingStepsizes.html>.
% * |seti.useTolIn|     : If it is set to 1 (default: 0), the inner 
%                         tolerance principle is used to stop the inner 
%                         iterations (primal-dual algorithm), see 
%                         <minPda.html> and <minTolIn.html>.
% * |seti.useTolOut|    : If it is set to 1 (default: 0), the outer 
%                         tolerance principle is used to stop the inner 
%                         iterations, see <minPda.html>, <minTolOut.html>.
%
% In case of |seti.inv = 'pda'| and |seti.useTolIn = 1| and the following
% fields does not exist, default values are set in this function.
%
% * |seti.ThetaStart|   : (default: 0.925), see <minTolIn.html>.
% * |seti.ThetaMax|     : (default: 0.95), see <minTolIn.html>.
% * |seti.TolGamma|     : (default: 0.90), see <minTolIn.html>.
% * |seti.pdaNmax|      : Maximal iteration number, if |seti.useTolIn| or
%                         |seti.useTolOut| are 1 (default : 250).
%
% In case of |seti.inv = 'pda'| and |seti.useTolOut = 1| and the following
% fields does not exist, default values are set in this function.
%
% * |seti.relDisTol|    : Outer tolerance (default: 0.05), see <minTolOut.html>.
% * |seti.pdaNmax|      : Maximal iteration number, if |seti.useTolIn| or
%                         |seti.useTolOut| are 1 (default : 250).
%
% In case of |seti.inv = 'pda'| and the usage of |fp2| (see <setInvType.html>),
% that is not supported in the public versions it is checked if the following 
% fields exist and are not empty. Otherwise the program is stopped.
%
% * |seti.obsMask|      : Mask of pixels in ROI the assumed obstacle 
%                         (e.g. by factorization method)
%                         stored as logical vector of size seti.nROI^seti.dim x 1.
%
% The following fields should have been set in <setInvType.html>:
%
% * |seti.inv|
% * |seti.alpha|
% * |seti.beta|
% * |seti.pNorm|
% * |seti.qNorm|
% * |seti.tau|
%
% See <recon.html> for the following fields of |seti|:
%
% * |seti.saveqROIcomp|
% * |seti.loadqROIcomp|
% * |seti.savedata|
%
% Currently unused fields:
%
% * |seti.lambda|   : Wavelength (default: $\lambda = 2\pi/k$ with 
%                     wave number $k = \texttt{seti.k}$) (currently not used anywhere).
%
%% Output Argument
%
% * |seti|  : structure array, see "Optional Input Arguments of structure
% |seti|", because default values are set, if fields was not already
% defined.
%
% * |seti.usefp2|       : 1 if |fp2| is used in |seti.tF| or |seti.tG| (otherwise 0)
%                         (not supported in public version.)
%
%% See Also
% * <recon.html>
% * <minPda.html>
% * <pda.html>
%
%% Code
%
function seti = checkConsisRec(seti,dispDepth)

%%
% *setWavelet*
%
% Wavelets are not supported in public version.

seti = checkfield(seti,'useWavelet',0,dispDepth);

if seti.useWavelet == 1

    seti = checkfield(seti,'wavWeight',0,dispDepth);
    seti = checkfield(seti,'genWhiteNoiseNew',0,dispDepth);

    if seti.dim == 2
        seti = checkfield(seti,'adaptWeightsMan',1,dispDepth);
    elseif seti.dim == 3
        seti = checkfield(seti,'adaptWeightsMan',0,dispDepth);
        disp('Note: Parameter "adaptWeightsMan" = 1 is not supported in case of seti.dim = 3 yet.')
    else
        error('Dimension has to be 2 or 3, see checkConsisRec.m')
    end

    seti = checkfield(seti,'wavelet','cdf97',dispDepth);

    seti = checkfield(seti,'extension','symmetric',dispDepth); % extension: symmetric or zeropadded
    
    if ~( isfield(seti,'smin') && ...
            (seti.smin == round(seti.smin)) && (seti.smin >= 8) )
        seti.smin = 8;
        setmessage(seti,'smin',dispDepth);
    end
   
    if ~( isfield(seti,'samples') && ...
            (seti.samples == round(seti.samples)) && (seti.samples >= 1000) )
        seti.samples = 10000;
        setmessage(seti,'samples',dispDepth);
    end
    
%%
% *Example for a complete but long consistency check in case of seti.wavelet*

%     if ~( isfield(seti,'wavelet') && ...
%             ( strcmp(seti.wavelet,'cdf97') || strcmp(seti.wavelet,'haar') || strcmp(seti.wavelet,'legall53') ) )
%         seti.wavelet = 'cdf97';
%         setmessage(seti,'wavelet');
%     end

end

%%
% *Variational reconstruction process*
%
% For reference see
%
% * <recon.html>
% * <setInvType.html>
% * <pda.html>
% * <setFuncsPda.html>

seti = checkfield(seti,'inv','pda',dispDepth);
seti = checkfield(seti,'recname','contrast',dispDepth); % recname: contrast or source

% -

seti = checkfield(seti,'alpha',500,dispDepth);
seti = checkfield(seti,'beta',1E-5,dispDepth);

seti = checkfield(seti,'pNorm',2,dispDepth);
seti = checkfield(seti,'qNorm',1,dispDepth);

seti = checkfield(seti,'physBounds',[-1,3,0,3],dispDepth);

if seti.useWavelet == 1 && strcmp(seti.inv,'shrinkage')
    seti.physBounds = [-inf,inf,-inf,inf]; % no bounds for wavelet coefficients
    % in case of pda physBounds can be used in wavelet case too
    disp('Wavelets are used and shrinkage is used...');
    disp('... so set parameter "physBounds" to [-inf,inf,-inf,inf].');
end

seti = checkfield(seti,'tolDelta',1E-14,dispDepth); % in case of pda
% tolDelta in setFuncsPda used in function M2 = G.

% -

seti = checkfield(seti,'tau',1.1,dispDepth);
seti = checkfield(seti,'useDis',1,dispDepth);

% -

seti = checkfield(seti,'nOut',30,dispDepth);

%%
% *Reconstruction with primal-dual algorithm*
%
% For reference see
%
% * <minPda.html>
% * <pda.html>

if strcmp(seti.inv,'pda')
    seti = checkfield(seti,'pdaN',50,dispDepth);
    seti = checkfield(seti,'pdaStepsize','fix',dispDepth); % pdaStepsize: fix or adaptive
    seti = checkfield(seti,'vartheta',0.5,dispDepth); % used in adaptive stepsize...

    % Check if the background indices are defined if seti.tF or seti.tG contains fp2
    if sum([strcmp('fp2',seti.tF), strcmp('fp2',seti.tG)]) == 1
        % fp2 is used in seti.tF or seti.tG
        % If seti.obsMask is not defined or empty, you should use fp instead of fp2.
        if ~isfield(seti,'obsMask') || isempty(seti.obsMask) || size(unique(seti.obsMask),1) ~= 2
            disp('The usage of fp2 (set background pixels to zero) instead of fp does not make sense,')
            disp('because one of the following reasons:')
            disp('1. seti.obsMask is not defined.')
            disp('2. seti.obsMask is empty.')
            disp('3. seti.obsMask is the whole ROI or nothing of ROI.')
            error('You should use fp, e. g. setting seti.invNo to 6.')
        end
        seti.usefp2 = 1;
    else
        seti.usefp2 = 0;
    end
    
    % tolerance outside pda:
    seti = checkfield(seti,'useTolOut',0,dispDepth);
    if seti.useTolOut == 1
        seti = checkfield(seti,'relDisTol',0.05,dispDepth);
        
        if ~isfield(seti,'pdaNmax')
            seti.pdaNmax = 250;
            setmessage(seti,'pdaNmax',dispDepth);
        end

        if isfield(seti,'pdaNmax')
            seti.pdaN = seti.pdaNmax;
            if dispDepth >= 1
                disp('Because useTolIn == 1, parameter pdaN is overwritten by "pdaNmax".');
            end
        end

    end

    % tolerance inside pda:
    seti = checkfield(seti,'useTolIn',0,dispDepth);
    if seti.useTolIn == 1
        seti = checkfield(seti,'ThetaStart',0.925,dispDepth);
        seti = checkfield(seti,'ThetaMax',0.95,dispDepth);
        seti = checkfield(seti,'TolGamma',0.90,dispDepth);

        if ~isfield(seti,'pdaNmax')
            seti.pdaNmax = 250;
            setmessage(seti,'pdaNmax',dispDepth);
        end

        if isfield(seti,'pdaNmax')
            seti.pdaN = seti.pdaNmax;
            if dispDepth >= 1
                disp('Because useTolIn == 1, parameter pdaN is overwritten by "pdaNmax".');
            end
        end
    end
end

%%
% *Reconstruction with shrinkage*
%
% Reconstruction with shrinkage is not supported in public version.

if strcmp(seti.inv,'shrinkage')
    seti = checkfield(seti,'stepsizeStart',1,dispDepth);
    seti = checkfield(seti,'useBarBor',0,dispDepth);
    seti = checkfield(seti,'useProjStepsize',0,dispDepth);
    seti = checkfield(seti,'useArmijo',0,dispDepth);
    seti = checkfield(seti,'maxArmijo',10,dispDepth);
end

%%
% *Save and load previously computed qROI*
%
% See <recon.html> for reference.

if isfield(seti,'dirname') && isfield(seti,'fileSuffix')
    seti = checkfield(seti,'saveqROIcomp',sprintf('%s/save_qROIcomp_iOutStop%s.mat',seti.dirname,seti.fileSuffix),dispDepth);
    %in qROIsaveFile is fileSuffix important in case of varalpha (various alpha) e.g.
else
    if dispDepth >= 1
        disp('   The field saveqROIcomp of structural array seti was not defined,')
        disp('   because the fields dirname or/and fileSuffix are not defined.')
    end
end

seti = checkfield(seti,'loadqROIcomp',[],dispDepth); % default: [], so no load...

%%
% *Save reconstructed contrast, ...*
%
% See <recon.html> for reference.

seti = checkfield(seti,'savedata',0,dispDepth);

%%
% *Currently unused...*
%
% Computes the wavelength $\lambda = 2\,\pi/k$ with wavenumber $k$.

seti = checkfield(seti,'lambda',2*pi/seti.k,dispDepth);

end
