%% setRecon
% Sets required fields in structure array |seti| for reconstruction
% process in <recon.html>.
%
%% Syntax
%
%   seti = setRecon(seti)
%   seti = setRecon(seti,dispDepth)
%
%% Input Arguments
%
% * seti    :   structure array
%
% *Optional Input Argument*
%
% * dispDepth   : Depth of displayed messages (0: no, 1 or greater: yes).
%
%% Output Arguments
%
% * seti    :   structure array
%
% The added fields in |seti| can be found in the called functions
%
% * <setInvType.html>
% * <checkConsisRec.html>
% * <shrinkFuncs.html>
% * <setFuncsPda.html>
%
%% More About
%
% * |shrinkFuncs| was originally written for reconstruction process by 
%   shrinkage, but the defined |seti.shrkRe| and |seti.shrkIm| are also 
%   needed in <pda.m>.
%
% *Notation in reconstruction process*
%
% * *Tikhonov functional* is called |MT|:
% * |MT = M1 + M2|.
% * MT: Tikhonov functional with parts M1 and M2.
% * M stands for minimization.
%
% * *Forward operator* is called |FF|:
% * $\texttt{FF} = \mathcal{F}$.
% * FF is the contrast-to-measurement operator.
% * Useful to call it $F_\mathrm{meas}$, 
%   because then clearly the measurements (receivers) are meant
%   (somtimes we evalute on ROI).
% * |FF(q) = FFqMeas (= uScattRX)| is on receivers (measurements):
% * ---- $\texttt{FmeasExact}$ is the exact scattered field at receivers positions.
% * ---- $\texttt{FmeasDelta} = F_\mathrm{meas}^\delta$ is the scattered 
% field with noise at receivers' positions.
% * (|FmeasComp| would be a consequently choosen name for currently
% computed, but is not used. Instead we use FFqMeas.)
% * (Analog would be consequently: |FROIcomp|, |FROIexact|.)
% * |FFqmF = FFqMeas - FmeasDelta|
% * Note that in |FFqmF| the character |m| stands for minus.
%
% * *Forward operator on ROI*:
% * FFqROI (= uScattROI) (analog to $\mathcal{F}(q)$, 
%   but on ROI and not not on receivers positions.)
%
% * *Frechet derivative of F* is called |DFF|
% * $\texttt{DFF} = \mathcal{F}'(q)[\cdot]$.
%
% * *Jacobi matrices* are called |JA| and |JB|:
% * $\texttt{DFFqh} = \mathcal{F}'(q)[h] = \texttt{JA*diag(h)*JB}$.
%
% * *Adjoint of derivative* is called |ADFF|:
% * $\texttt{ADFFq} =
%   [\mathcal{F}'(q)]^\ast[\mathcal{F}(q)-F_\mathrm{meas}^\delta]$
%   $\quad\quad$ (in case of no wavelet)
% * $\texttt{ADFFq} = \texttt{iWstar}
%   [\mathcal{F}'(q)]^*[\mathcal{F}(q)-F_\mathrm{meas}^\delta]$
%   $\quad\quad$ (in case of wavelet, i.e. useWavelet = 1, which is not
%   available in public version.)
%
%% See Also
% * <setInvType.html>
% * <checkConsisRec.html>
% * <shrinkFuncs.html>
% * <setFuncsPda.html>
%
function seti = setRecon(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

%% sets for reconstruction process

if dispDepth >= 1
    disp('-- setInvType --')
end
    seti = setInvType(seti,dispDepth);

if dispDepth >= 1
    disp('-- check consistency of parameters used in inversion process --')
end
seti = checkConsisRec(seti,dispDepth);

seti = checkfield(seti,'closed',0); % is seti.closed was not set in setInput set it to 0.
if seti.closed == 1
    if dispDepth >= 1
        disp('-- setWavelet --')
    end
    seti = setWavelet(seti);
end

%% Definition of functions

if dispDepth >= 1
    disp('-- set shrinkage functions --')
end
seti = shrinkFuncs(seti); % defines seti.shrkRe and seti.shrkIm that are also needed in pda

if strcmp(seti.inv,'shrinkage')
    if dispDepth >= 1
        disp('-- setFuncsShrink --')
    end
    seti = setFuncsShrink(seti); % not available in public version
elseif strcmp(seti.inv,'pda')
    if dispDepth >= 1
        disp('-- setFuncsPda --')
    end
    seti = setFuncsPda(seti);
end

%%
if isfield(seti,'qCDexact')
    seti = rmfield(seti,'qCDexact');
end

end

