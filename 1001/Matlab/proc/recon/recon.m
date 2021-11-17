%% recon
% Variational reconstruction of the contrast.
%
%% Syntax
%
%   seti = recon(seti)
%   seti = recon(seti,dispDepth)
%   seti = recon(seti,dispDepth,out)
%
%% Description
% |seti = recon(seti)| 
% does the variational reconstruction of the contrast by
% minimizing the Tikhonov functional 
% by _soft-shrinkage_ or _primal-dual algorithm_.
%
% |seti = recon(seti,dispDepth)| does the same as |seti =
% recon(seti)|, but the depth of displayed messages is controlled by
% |dispDepth|.
%
% |seti = recon(seti,dispDepth,out)| does the same as |seti = recon(seti,dispDepth)| in 
% case of |out = 0|, but plots the figures in case of |out = 1| and 
% additionally saves them in case of |out = 2|.
%
% * In public version only primal-dual algorithm is available.
% * The result is the contrast |seti.qROIcomp| after |seti.iOutStop| 
% outer iterations.
% * The reconstruction is plotted in case of |out = 1| and in case of 
%   |out = 2| additionally stored in the folder output.
%
% * The reconstruction is stopped after maximal |seti.nOut| iterations.
% * If |seti.useDis = 1| and a tolerance parameter |seti.tau > 1| was set,
% the reconstruction is stopped by *discrepancy principle*, i.e.
%
% stop if discrepancy $\leq$ $\texttt{seti.tau*seti.delta}$.
%
% * If |seti.loadqROIcomp| is not empty, 
% a path to a mat-file is expected, that contains 
% the reconstruction result |qROIcomp| after |iOutStop| outer iterations 
% from a previous computation. Then this file is loaded. 
% Further, the contrast and the stop index of outer iteration are set as 
% initial contrast |seti.qROIcomp| and iteration number |seti.iOutIni| 
% for current computation.
%
% * If |seti.saveqROIcomp| is not empty (i.e. |''| or |[]|), a filename to save 
% |qROIcomp| and |iOut| is expected. Default in <checkConsisRec.html> is
%
%  sprintf('%s/save_qROIcomp_iOutStop%s.mat',seti.dirname,seti.fileSuffix)
%
% (i.e. the mat-file is stored in the same folder as the figures).
% Then |seti.qROIcomp| and |seti.iOutStop| are saved as |qROIcomp| and
% |iOutStop|.
%
% * If |seti.savedata| is 1 (default: 0), 
%   relative discrepancies, errors, and differences of iterated contrasts
%   are saved in mat-files in the folder of the
%   figures (filenames are essentially save_dis.mat, save_err.mat, save_dif.mat.)
%
%
% *Structure of this file*
%
% # Optional: loads reconstruction result from previously computation.
% # Initialization.
% # Evaluation of previous outer iteration.
%   (Stop if maximal number of iterations is reached or 
%   discrepancy principle is fullfilled (optional)).
% # Plot results from previous outer iteration 
%   and stores them in folder output (optional, depends on argument |out|).
% # Compute next step (minimization by soft-shrinkage or primal-dual algorithm). Go to 3.
% # Optional: save discrepancy, error and reconstruction results.
%
%
%% Example
%
% *Example 1*
%
%   init;                   % Initialization (addpath...)
%   seti = struct;
%   seti = setData(seti);   % Set data: (experimental data, geometry, exact and noisy data)
%   seti = setRecon(seti);  % Settings for variational reconstruction
%   seti = recon(seti);     % Variational reconstruction (process)
%
% *Example 2*
%
% This example does the same as the first one, but displays more messages.
%
%   init;                     % Initialization (addpath...)
%   seti = struct;
%   seti = setData(seti,4);   % Set data: (experimental data, geometry, exact and noisy data)
%   seti = setRecon(seti,4);  % Settings for variational reconstruction
%   seti = recon(seti,4);     % Variational reconstruction (process)
%
% *Example 3*
%
% The example is essentially the code in start.m (<start.html>) 
% because we need previously definitions.
%
% Warning: variables are cleared and figures as well as files are saved
% without demand.
%
% The results are stored in the folder output.
%
%   init;                       % Initialization (addpath...)
%   setInput;                   % Set input and make directories for output
%   seti = setData(seti,4,2);   % Set data: (experimental data, geometry, exact and noisy data)
%   seti = setRecon(seti,4);    % Settings for variational reconstruction
%   seti = recon(seti,4,2);     % Variational reconstruction (process)
%
%% Input Arguments
%
% * |seti|    :   structure array
%
% Because this is a internal function we do not explain all fields in seti.
%
% * |seti.loadqROIcomp| : Path to load a mat-file containing the reconstructed 
%                         contrast |qROIcomp| after |iOutStop| outer iterations 
%                         (default: empty, i.\,e. |''| or |[]|, then no data is loaded).
% * |seti.saveqROIcomp| : Optional: Path to save reconstructed contrast
%                         (default: same folder as figures, i.e. 
%                         |sprintf('%s/save_qROIcomp_iOutStop%s.mat',seti.dirname,seti.fileSuffix)|)
%
% * |seti.nOut|         : Maximal number of outer iterations.
%
% * |seti.useDis|       : Use discrepancy principle to stop outer iteration?
%                         (If yes set it to 1.)
% * |seti.tau|          : Tolerance parameter for discrepancy principle
%                         (seti.tau > 1).
%
% * |seti.plotFreq|     : Frequency to plot figures of outer iteration steps
%                         (0: no figures).
% * |seti.savedata|     : If it is set to 1,
%                         relative discrepancies, errors, and differences 
%                         of iterated contrasts are saved
%                         (as |save_dis.mat|, |save_err.mat| and |save_dif.mat|).
%
% * |seti.qROIexact|    : Predefined (exact) contrast
%                         (complex vector of |size seti.nROI|^|seti.dim| x 1).
%
% *Further specific input parameters* for primal-dual algorithm can be found in
% <minPda.html> and <pda.html>.
%
% *Optional Input Argument*
%
% * |dispDepth|     : depth of displayed text.
%                     (from 0 to 5, default sparse: 2; default details: 4; too much information: 5);
% * |out|           : output depth for figures and files: no figure (0), plot figure (1), 
%                     save figures and files (2).
%                     (This argument requires that argument |dispDepth| was
%                     defined.)
%
%% Output Arguments
%
% * |seti|    :   structure array
%
% Because this is a internal function we do not explain all fields in seti.
%
% * |seti.qROIcomp|     : Reconstructed contrast
%                         (complex vector of size |seti.nROI|^|seti.dim| x 1).
% * |seti.iOutStop|     : Stop index of outer iterations.
%
%
% *Discrepancy, error, and difference*
%
% * |seti.dis|        : Relative discrepancy for each outer iteration
%                       (vector of size 1 x seti.nOut).
% * |seti.err|        : Relative error for each outer iteration
%                       (vector of size 1 x seti.nOut).
% * |seti.dif|        : Relative difference of iterated contrast qROI to previously qROI
%                       for each outer iteration
%                       (vector of size 1 x seti.nOut).
%
% Note that |seti.dis(iOut)| is the relative discrepancy 
% _after_ iteration |iOut|. (Analog seti.err(iOut) and seti.dif(iOut).)
%
% *Minimized Tikhonov functional*
%
% * |seti.MTv|        :   Result of minimized Tikhonov functional 
%                         for each outer iteration 
%                         (vector of size 1 x seti.nOut).
% * |seti.M1v|        :   Result of first part of min. functional
%                         for each outer iteration 
%                         (vector of size 1 x seti.nOut).
% * |seti.M2v|        :   Result of second part of min. functional 
%                         for each outer iteration 
%                         (vector of size 1 x seti.nOut).
%
% _More about M1v and M2v_
%
% * In case of primal-dual algorithm (i.e. |seti.inv = 'pda'|)
%   |M1v| is F and |M2| is G, see <pda.html> for F and G.
% * In case of shrinkage (i.e. |seti.inv = 'shrinkage'|)
%   (not available in public version) 
%   |M1| = discrepancy and |M2| = sparsity (setFuncsShrink.m)
%
% *structure array pdas*
%
% * |pdas|        : Structure array
%                   (only in case of pda, i.e. |seti.inv = 'pda'|).
%
% * |pdas.disLin| : Last relative discrepancy of linearized problem in inner iteration of pda
%                   for each outer iteration 
%                   (vector of size 1 x seti.nOut).
% * |pdas.relDis| : Quotient disLin/dis 
%                   for each outer iteration
%                   (rel. discrepancy of linearized problem / rel. discrepancy of the non-linearized problem)
%                   (vector of size 1 x seti.nOut).
% * |pdas.pdaNv|  : Number of inner iterations for each outer iteration
%                   (vector of size 1 x seti.nOut).
% * |pdas.ThetaiOutV| : Inner tolerance for each outer iteration
%                       (in case of inner tolerance principle, see
%                       <minTolIn.html>).
%
% *Further output arguments*
%
% * |seti.disIni|   : Relative discrepancy before 1st iteration.
% * |seti.errIni|   : Relative error before 1st iteration.
%
%
%% Output Figures
%
% The functions in the code, that plot figures are:
%
% * subplots            :   figure 11
% * plotAndSaveFigures  :   figures 12-16 and in 3D case figures 21-26
% * pdaPlot             :   figures 31--36
%
% A description of all figures in the package are in <start.html>, Section
% "Output: Figures".
%
% Note that the creation and saving of plots also depends on the input
% argument |out|.
%
%
%% More About
%
% The outer iteration of the variational reconstruction process is
% briefly described in <start.html>, Section "More About". More information
% are available in [1].
%
% The primal-dual algorithm and its application to the inverse scattering
% problem is briefly described in <pda.html>, Section "More About". For
% further information of application see Section 4 in [1]. For details in
% context of primal-dual algorithm see [2].
%
% The initial contrast |seti.qROIcomp| is set to zero if no contrast is
% loaded from a previous reconstruction. 
%
%
%% References
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
% * [2] Antonin Chambolle and Thomas Pock. A first-order primal-dual algorithm for convex problems with applications to imaging. _Journal of Mathematical Imaging and Vision_, 40(1):120-145, 2011.
%
%% See Also
%
% * <start.html>
% * <pda.html>
% * <setFuncsPda.html>
% * <minPda.html>
% * <minTolIn.html>
%
%% Code
%
function seti = recon(seti,varargin)

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
% *Optional: loads reconstruction result from a previous computation.*
%
if ~isempty(seti.loadqROIcomp)
    disp('# Loading qROI from previous calculation.')
% -- load old: start --
%     load(sprintf('%s',seti.loadqROIcomp)); % qROIcomp and iOutStop are inside
%     seti.qROIcomp = qROIcomp;
%     seti.iOutIni = iOutStop;
%     clear qROIcomp iOutStop;
% -- load old: end --
% -- load new: start --
    sloaded = load(sprintf('%s',seti.loadqROIcomp)); % qROIcomp and iOutStop are inside (sloaded: struct)
    seti.qROIcomp = sloaded.qROIcomp;
    seti.iOutIni = sloaded.iOutStop;
    clear sloaded;
% -- load new: end --
else
    seti.qROIcomp = zeros(size(seti.qROIexact));
    seti.iOutIni = 0;
end

% FmeasDelta = seti.FmeasDelta; % perturbed data.
% If you call pertuped data, use consequently seti.FmeasDelta,
% if you are not in a function.

%%

% -- iteration - start with q_0 = 0 or loaded qROI
seti = minimization(seti,out,dispDepth);

end

%% Code: subfunction: minimization
%
function seti = minimization(seti,out,dispDepth)

if dispDepth >= 1
    disp(' ')
    fprintf('   Remind: alpha = %g\n',seti.alpha);
    fprintf('   Remind: beta  = %g\n',seti.beta);
end

% qROIexact = seti.qROIexact; % use consequently seti.qROIexact
% qROIcomp = seti.qROIcomp; % use consequently seti.qROIcomp

seti.dis = zeros(1,seti.nOut); % seti.dis(iOut) = relative discrepancy after iOut
seti.err = zeros(1,seti.nOut); % seti.err(iOut) = relative error after iOut
seti.dif = zeros(1,seti.nOut); % seti.dif(iOut) = relative difference of iterated qROI to previously qROI after iOut

seti.MTv = zeros(1,seti.nOut); % results of minimized Tikhonov functional
seti.M1v = zeros(1,seti.nOut); % results of first part of min. functional
seti.M2v = zeros(1,seti.nOut); % results of second part of min. functional

if strcmp(seti.inv,'pda')
    pdas.disLin = ones(1,seti.nOut); % store results of last relative lin. discrepancy in inner iteration of pda
    pdas.relDis = zeros(1,seti.nOut); % disLin/dis (not in pda, but outside)

    pdas.pdaNv = zeros(1,seti.nOut); % number of inner iterations in each outer iteration
    pdas.ThetaiOutV = zeros(1,seti.nOut);
end

%%
% *minimization: Initialization*
%
ticMin = tic;

iOut = seti.iOutIni; % outer iteration
% (default: 0; higher if previously computed qROI is used; is set in start.m)

% initial input is seti.qROIcomp
if dispDepth >= 1
    disp(' ')
    fprintf(' - iOut = %02d (initial, no reconstruction done)\n',iOut); % initial...
end
if iOut ~= 0
    FFqMeas = mimo(seti, seti.qROIcomp, 'simo'); % only needed to compute dis
    seti.disIni = compDis(seti,FFqMeas,seti.FmeasDelta); % relative discrepancy before 1st iteration
    seti.errIni = compErr(seti.qROIexact,seti.qROIcomp); % relative error before 1st iteration
else
    if dispDepth >= 1
        disp('   Because iOut == 0 set disIni = 1 and errIni = 1 (otherwise compute them).')
    end
    seti.disIni = 1;
    seti.errIni = 1;
end
if dispDepth >= 1
    fprintf('   disIni = %g \n',seti.disIni)
    fprintf('   errIni = %g \n',seti.errIni)
end

if ~exist('ilast','var')
    ilast = false;
end

while 1
    %% 
    % *minimization: Evaluation of previous outer iteration*

    % qROIcomp: real and imag
    if dispDepth >= 1
        fprintf('   qROIcomp: real [%+3.2f,%+3.2f], imag [%+3.2f,%+3.2f].\n',...
            min(real(seti.qROIcomp)),max(real(seti.qROIcomp)),min(imag(seti.qROIcomp)),max(imag(seti.qROIcomp)));
        fprintf('   qROIexact: real [%+3.2f,%+3.2f], imag [%+3.2f,%+3.2f].\n',...
            min(real(seti.qROIexact)),max(real(seti.qROIexact)),min(imag(seti.qROIexact)),max(imag(seti.qROIexact)));
    end
    
    if iOut == seti.iOutIni
        dis = seti.disIni;
    else
        dis = seti.dis(iOut);
    end
    [ibreak,seti.iOutStop] = reconStop(seti,iOut,dis,dispDepth);

    if ilast == true
        seti.iOutStop = iOut;
        ibreak = true;
    elseif ilast == false && isfield(seti,'pdaNlast') && seti.pdaNlast > 0
        if iOut == seti.nOut-1 % set ilast = true because last iteration will follow
            disp(' ')
            disp('   seti.pdaNlast > 0. The last iteration will done by a higher number of inner iterations.')
            seti.pdaN = seti.pdaNlast; % pdaN is overwritten for the last step.
            ibreak = false;
            ilast = true;
        elseif ibreak == true % in the case of discrepancy principle
            disp(' ');
            fprintf('   Start the last outer iteration with %g inner iterations.\n', seti.pdaNlast)
            seti.pdaN = seti.pdaNlast; % pdaN is overwritten for the last step.
            ibreak = false;
            ilast = true;
        end
    end


    %% 
    % *minimization: Plot results from previous outer iteration*
    %
    % The figures are described in <start.html>.
    %
    if (out >= 1 && seti.plotFreq ~= 0 &&...
            (iOut == 1 || ibreak == true || (isfield(seti,'pdaNlast') && seti.pdaNlast > 0 && iOut == seti.nOut-1)...
             || ilast == true || floor(iOut/seti.plotFreq) == iOut/seti.plotFreq))
    % condition "iOut == seti.nOut-1" etc. to plot the step before the last step (interesting in the case of pdaNlast...) 
    
        if iOut ~= seti.iOutIni
            subplots(seti,seti.qROIexact,seti.qROIcomp,iOut); % figure 11
        end
        plotAndSaveFigures(seti,seti.qROIexact,seti.qROIcomp,iOut,out); % save figure 11; generate figures 12--16 and in 3D case 21--26
        
        if iOut > 0 && strcmp(seti.inv,'pda') && isfield(pdas,'pdaStopInd')
            pdaPlot(iOut,seti,pdas,out) % figures 31--36
        end
    end

    % time
    if iOut ~= seti.iOutIni
        tocOut = toc(ticOut);
        if dispDepth >= 1
            disp(' ')
            fprintf('   Elapsed time of outer minimization iOut = %03d is %05.1f min.\n',iOut,tocOut/60);
        end
        
        tocMin = toc(ticMin);
        if dispDepth >= 1
            disp(' ')
            fprintf('   Elapsed time of overall minimization process is %05.1f min.\n',tocMin/60);
        end
    end
    
    if ibreak == true
        break;
    end
    % break does mean: result is seti.qROIcomp after iOut iterations

    %% 
    % *minimization: Compute next step*
    %
    % (The new contrast "qROInew" is seti.qROIcomp after minShrink or minPda.)

    ticOut = tic;
    iOut = iOut + 1;
    if dispDepth >= 1
        disp(' ')
        fprintf(' - iOut = %02d\n',iOut);
    end
    
    switch seti.inv
        case 'shrinkage'
            % Note that the case seti.inv = 'shrinkage' is not available in
            % public version.
            if iOut == seti.iOutIni+1
                ADFFqPrev = zeros(seti.nROI^seti.dim,1); % to have a value
                qROIcompPrev = zeros(seti.nROI^seti.dim,1); % to have a value
            end
            % ADFFqPrev is needed in minShrink in case of Barzilai-Borwein
            [seti,FFqMeas,qROIcompPrev,ADFFqPrev] = minShrink(seti,iOut,seti.qROIcomp,qROIcompPrev,ADFFqPrev);
        case 'pda'
            if dispDepth >= 2
                disp('    - Start minimization by PDA')
            end
            if iOut == seti.iOutIni + 1 % otherwise it is done after pda algorithm
                FFqMeas = mimo(seti, seti.qROIcomp,'simo'); % needs upscaled qROI
            end
            qROIcompPrev = seti.qROIcomp;
            [seti,FFqMeas,pdas] = minPda(seti,iOut,seti.qROIcomp,pdas,FFqMeas,dispDepth);
    end
    % output: seti.qROIcomp (new reconstructed contrast after step iOut)
    
    % *Compute discrepancy, error and difference*
    %
    seti = disErrDif(FFqMeas,qROIcompPrev,seti,iOut,dispDepth);

    % *Compute quotient disLin/dis*
    % disLin/dis = rel. discrepancy of linearized problem / rel. discrepancy of the non-linearized problem
    %
    if iOut > 0 && strcmp(seti.inv,'pda')
        pdas.relDis(iOut) = pdas.disLin(iOut)/seti.dis(iOut);
        if dispDepth >= 3
            fprintf('   relDis = %5.3f \n',pdas.relDis(iOut));
        end
    end

end

%%
% *minimization: Optional: save relative discrepancies, errors, and differences of iterated contrasts*
%
if out == 2
    savedata(seti,seti.dis,seti.err,seti.dif);
end

% *Optional: save reconstructed contrast qROIcomp and stop index iOutStop
% for further computations*
%
if out == 2 && isfield(seti,'saveqROIcomp') && ~isempty(seti.saveqROIcomp)
    saveqROI(seti);
end

end

%% Code: Further Subfunctions

%%
% *subfunction: savedata*
%
% Save relative discrepancies, errors, and differences of iterated contrasts.
%
function savedata(seti,dis,err,dif)

if seti.savedata
    % Remove zeros from dis, err, and dif.
    % (dis ~= 0) or 1:iOutStop... should work too...
    disNonZeros = dis(dis ~= 0); %#ok (suppress MATLAB warning, because will be saved.)
    errNonZeros = err(err ~= 0); %#ok (suppress MATLAB warning, because will be saved.)
    difNonZeros = dif(dif ~= 0); %#ok (suppress MATLAB warning, because will be saved.)
    %  Save relative discrepancies, errors, and differences of iterated contrasts.
    save(sprintf('%s/save_dis%s.mat',seti.dirname,seti.fileSuffix),'disNonZeros');
    save(sprintf('%s/save_err%s.mat',seti.dirname,seti.fileSuffix),'errNonZeros');
    save(sprintf('%s/save_dif%s.mat',seti.dirname,seti.fileSuffix),'difNonZeros');
    clear disNonZeros;
    clear errNonZeros;
    clear difNonZeros;
end
end

%%
% *Subfunction: reconStop*
%
% Stops the reconstruction process, when maximal number of outer iterations
% |seti.nOut| is reached or discrepancy principle is fullfilled (optional).
%
function [ibreak,iOutStop] = reconStop(seti,iOut,dis,dispDepth)

if iOut == seti.nOut
    iOutStop = iOut;
    if dispDepth >= 1
        disp(' ')
        fprintf('   Reconstruction stopped after %0d outer iterations (max. number reached).\n',iOutStop)
    end
    ibreak = true;
end

if iOut > 0 && isnan(dis)
    iOutStop = iOut;
    if dispDepth >= 1
        disp(' ')
        fprintf('   Reconstruction stopped after %0d outer iterations (discrepancy not a number).\n',iOutStop);
    end
    ibreak = true;
end

if iOut > 0 && seti.useDis == 1 && (dis <= seti.tau*seti.delta)
    iOutStop = iOut;
    if dispDepth >= 1
        disp(' ')
        fprintf('   Reconstruction stopped after %0d outer iterations (by discrepancy principle).\n',iOutStop);
    end
    ibreak = true;
end

if ~exist('ibreak','var')
    ibreak = false;
end

if ~exist('iOutStop','var')
    iOutStop = [];
end

end

%%
% *Subfunction: compDis*
%
% Computes relative discrepancy
%
function dis = compDis(seti,FFqMeas,FmeasDelta)
dis = normws(FFqMeas-FmeasDelta,seti)/normws(FmeasDelta,seti); % relative discrepancy
end

%%
% *Subfunction: compErr*
%
% Computes relative error
%
function err = compErr(qROIexact,qROIcomp)
err = norm(qROIexact-qROIcomp,2)/norm(qROIexact,2); % relative error in 2-Norm
end

%%
% *Subfunction: disErrDif*
%
% Computes, stores and displays:
%
% * realtive discrepancy,
% * relative error, and
% * relative difference of iterated contrast to previous one
%

function seti = disErrDif(FFqMeas,qROIcompPrev,seti,iOut,dispDepth)
% relative discrepancy
seti.dis(iOut) = compDis(seti,FFqMeas,seti.FmeasDelta);
clear FFqmF;
% relative error
seti.err(iOut) = compErr(seti.qROIexact,seti.qROIcomp);
% relative difference of iterated qROI to previously qROI
seti.dif(iOut) = normroi(seti.qROIcomp-qROIcompPrev,seti)/normroi(seti.qROIcomp,seti);
% display:
if dispDepth >= 2
    fprintf('   dis = %5.3f \n',seti.dis(iOut))
    fprintf('   err = %5.3f \n',seti.err(iOut))
    fprintf('   dif = %5.3f \n',seti.dif(iOut))
end
end

%%
% *Subfunction: saveqROI*
%
% Saves the reconstructed contrast |qROIcomp| and the stop index |iOutStop|
% in a file |seti.saveqROIcomp|.
%
function seti = saveqROI(seti)
    disp(' ')
    disp('# Save qROIcomp for further calculations')
    qROIcomp = seti.qROIcomp; %#ok (suppress MATLAB warning, because will be saved.)
    iOutStop = seti.iOutStop; %#ok (suppress MATLAB warning, because will be saved.)
    save(seti.saveqROIcomp,'qROIcomp','iOutStop'); % save qROIcomp and iOutStop
    clear qROIcomp iOutStop;
end
