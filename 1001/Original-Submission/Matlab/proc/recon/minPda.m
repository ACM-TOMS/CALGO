%% minPda
%
% This function is called in recon.m and minimizes the non-linear problem
% by linearizing it and solving the linearized problem by minimization with
% primal-dual algorithm in pda.m.
%
%% Syntax
%
%   [seti,FFqMeas,pdas] = minPda(seti,iOut,qROIcomp,pdas,FFqMeas,dispDepth)
%
%% Description
%
% |[seti,FFqMeas,pdas] = minPda(seti,iOut,qROIcomp,pdas,FFqMeas,dispDepth)|
% computes the outer iteration step |iOut| 
% to improve the reconstruction of the contrast |qROIcomp|, 
% that results in measurement data |FFqMeas|,
% from given data |seti.FmeasDelta|.
% With |dispDepth| the frequency of displayed text messages is controlled.
% 
% The reconstruction is done by linearizing the non-linear problem and
% solve the linearized one by minimization with primal-dual algorithm.
%
% The resulting contrast is saved in |seti.qROIcomp| and the corresponding
% measurement data in |FFqMeas|.
% 
% * If |seti.useTolIn = 1| the inner tolerance principle is used to stop the
% inner iterations (primal-dual algorithm).
% * If |seti.useTolOut = 1| the outer tolerance principle is used to stop the
% inner iterations.
% * The primal-dual algorithm stops at the latest afer |seti.pdaN|
% iteration steps.
%
%
%% Input Arguments
%
% Specific input arguments for primal-dual algorithm are described in
% <pda.html>.
%
% * seti            :   structure array
% * seti.pdaNmax    :   maximal iteration number of inner iteration,
%                       if useTolIn or useTolOut (default : 250)
%
% * iOut        :   iteration index of outer iteration
% * qROIcomp    :   reconstruction of the contrast
%                   (complex vector of size seti.nROI^seti.dim x 1)
% * pdas        :   structure array for pda specific outputs, 
%                   see below and also <recon.html>.
% * |FFqMeas|  :    scattered field evaluated on receivers positions
%                   (complex matrix of size seti.measNb x seti.incNb)
%
% * |dispDepth| : depth of displayed test messages (greater or equal 2 to
% get messages).
%
% *structure array pdas*
%
% * pdas.pdaStopInd     :   stop index of inner iteration (i.e. primal-dual algorithm)
% * pdas.MTvN           :   result of minimized Tikhonov functional 
% * pdas.M1vN           :   result of first part of min. functional
% * pdas.M2vN           :   result of second part of min. functional
% * See also <recon.html>.
%
% * pdas.relLinDisInPda :   quotient: relLinDis = disLin / dis (vector of size seti.pdaN x 1)
% * pdas.disLinInPda    :   discrepancy of linearized problem for each inner iteration step (vector of size seti.pdaN x 1)
% * pdas.errInPda       :   relative error of the reconstructed contrast qROI (vector of size seti.pdaN x 1)
% * pdas.minf = minf    :   struct with parts of the minimization functional (Tikhonov functional)
% * See also <pda.html>.
%
% * pdas.disLin :   last relative discrepancy of linearized problem in inner iteration of pda
%                   for each outer iteration 
%                   (vector of size 1 x seti.nOut)
% * pdas.relDis :   quotient disLin/dis 
%                   for each outer iteration
%                   (rel. discrepancy of linearized problem / rel. discrepancy of the non-linearized problem)
%                   (vector of size 1 x seti.nOut)
% * pdas.pdaNv  :   number of inner iterations for each outer iteration
%                   (vector of size 1 x seti.nOut)
% * pdas.ThetaiOutV     :   inner tolerance for each outer iteration
%                           (in case of inner tolerance principle, see
%                           <minTolIn.html>.)
%
%
%% Output Arguments
%
% * seti    :   structure array
% * seti.qROIcomp   :   new computed reconstruction of the contrast
%                       (complex vector of size seti.nROI^seti.dim x 1)
% * FFqMeas :   scattered field of currently computed contrast evaluated on receivers positions
%               (complex matrix of size seti.measNb x seti.incNb)
% * pdas    :   structure array for pda specific output, 
%               expanded for currently computed iteration step 
%               (see "Input Arguments" for a reference).
%
%% More About
%
% # Inner iteration: pda minimizes the linearized Tikhonov functional, output h.
% # Update in outer iteration: $q := q + h$.
% # Compute scattered field of reconstructed contrast at receivers
% positions: $\mathcal{F}(q)$.
%
% See also <start.html> and [1] Section 4.
%
%% References
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <start.html>
% * <pda.html>
% * <minTolIn.html>
% * <minTolOut.html>
%
%% Code
%
function [seti,FFqMeas,pdas] = minPda(seti,iOut,qROIcomp,pdas,FFqMeas,dispDepth)
%%
% *Parts of minimization function*
%
% seti.M1 = F(K(h)), seti.M2 = G(h)

if iOut == 1 % to have a value...
    seti.M1v(iOut) = 0;
    seti.M2v(iOut) = 0;
    seti.MTv(iOut) = 0;
end

%%
% *Tolerance principle in pda (in [1] "Stopping strategy 2")*
%
if seti.useTolIn
    if iOut == seti.iOutIni + 1
        ThetaiOut = seti.ThetaStart;
    else
        ThetaiOut = pdas.ThetaiOutV(iOut-1);
    end
    pdas.ThetaiOutV(iOut) = minTolIn(ThetaiOut,pdas.pdaNv,iOut,FFqMeas,seti,dispDepth);
    clear ThetaiOut;
else
    pdas.ThetaiOutV(iOut) = 0; % to set a value...
end

%%
% *Tolerance principle outside of pda (in [1] "Stopping strategy 1")*
%
if seti.useTolOut
    if iOut > seti.iOutIni+1
        seti.pdaN = minTolOut(iOut,pdas.pdaNv,seti.dis,pdas.disLin,pdas.relDis,seti,dispDepth);
        % choice of pdaN for next outer iteration
    else
        seti.pdaN = seti.pdaN;
    end
end

%%
% *Minimization of linearized Tikhonov functional by PDA*
%
% * |qROI| is |qCVU|, i.e. complex vector upscaled
% * In pda we require *real* vector spaces. 
%   Only the output h is complex.
%
if dispDepth >= 2
    disp('    - Minimization of linearized Tikhonov functional by PDA (start)')
end
[h,pdaStopInd,MTvN,M1vN,M2vN,relLinDisInPda,disLinInPda,errInPda,minf] = pda(iOut,qROIcomp,pdas.ThetaiOutV(iOut),seti,dispDepth); %FCGC only to plot
% Remember in pda: DFFq = @(xnRVD) JA*diag(seti.GU(seti.T(xnRVD)))*JB;
if dispDepth >= 2
    disp('    - Minimization of linearized Tikhonov functional by PDA (end)')
end

%%
% *Store output of pda in struct pda (pda struct: pdas) (except h)*
%
pdas.pdaStopInd = pdaStopInd;
pdas.MTvN = MTvN;
pdas.M1vN = M1vN;
pdas.M2vN = M2vN;
pdas.relLinDisInPda = relLinDisInPda;
pdas.disLinInPda = disLinInPda;
pdas.errInPda = errInPda;
pdas.minf = minf;

clear pdaStopInd MTvN M1vN M2vN relLinDisInPda disLinInPda errInPda minf;

%%
% *Store further values in struct pdas and seti*
%

% COMMENTED OUT: sigmaVal...
%
%size(sigmaVal) e.g. 25 x 1
%pdaStopInd = length(sigmaVal)
% can happen that sigmaVal is too small for sigmaMat
% fill the rest with zeros.
% does not work in case of useTolOut because pdaN changes...:
% if pdaStopInd < seti.pdaN
%    sigmaVal(seti.pdaN) = 0; % last is set 0, so the rest is filled with zero
%    tauVal(seti.pdaN) = 0;
%end

pdas.disLin(iOut) = pdas.disLinInPda(pdas.pdaStopInd); % lin. discrepancy

%pdaNv(iOut) = seti.pdaN; % does not respect inner tolerance principle
pdas.pdaNv(iOut) = pdas.pdaStopInd;
pdas.pdaStopInd = pdas.pdaStopInd;

%...N: vector of length N (from internal iteration)
seti.MTv(iOut) = pdas.MTvN(pdas.pdaStopInd);
seti.M1v(iOut) = pdas.M1vN(pdas.pdaStopInd);
seti.M2v(iOut) = pdas.M2vN(pdas.pdaStopInd);
if dispDepth >= 3
    fprintf('   MTv = %g, M1v = %g, M2v = %g\n',seti.MTv(iOut),seti.M1v(iOut),seti.M2v(iOut));
end

%%
% *Update in outer iteration*
if dispDepth >= 2
    disp('    - Update in outer iteration')
end
seti.qROIcomp = qROIcomp + h; % complex vector upscaled

%%
% *Compute scattered field of reconstructed contrast at receivers positions*
if dispDepth >= 2
    disp('    - Compute FF(q) for discrepancy and next step')
end
FFqMeas = mimo(seti, seti.qROIcomp, 'simo'); % needs upscaled qROI

end
