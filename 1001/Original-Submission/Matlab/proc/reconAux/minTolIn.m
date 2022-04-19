%% minTolIn
% Inner tolerance principle to stop the primal-dual algorithm.
% 
%% Syntax
%
%   ThetaiOut = minTolIn(ThetaiOut,pdaNv,iOut,FFq,seti)
%   ThetaiOut = minTolIn(ThetaiOut,pdaNv,iOut,FFq,seti,dispDepth)
%
%% Description
% |ThetaiOut = minTolIn(ThetaiOut,pdaNv,iOut,FFq,seti)|
% computes the inner tolerance |ThetaiOut|. See "More About" for details.
% |ThetaiOut = minTolIn(ThetaiOut,pdaNv,iOut,FFq,seti,dispDepth)| does the same but allows to control
% the displayed messages by |dispDepth|.
%
% To use it in context of computational framework set seti.useTolIn = 1.
%
%% Example
%
% See <minPda.html> to see the usage in computational framework.
%
% Here we give a small example.
%
%   init;
%
%   seti.ThetaStart = 0.925
%   seti.ThetaMax = 0.95
%   seti.TolGamma = 0.90
%
%   ThetaiOut = seti.ThetaStart;
%
%   seti.nOut = 30; % maximal number of outer iterations
%   pdaNv = zeros(1,seti.nOut);
% 
%   iOut = 1;
%
%   seti.incNb = 2;  % number of transmitters
%   seti.measNb = 5; % number of receivers
%   
%   FFq = rand(seti.measNb, seti.incNb);
%
%   seti.tau = 2.5;
%   seti.delta = 0.01;
%   seti.FmeasDelta = rand(seti.measNb,seti.incNb);
%
%   seti.pNorm = 2; % pNorm..., see normws.m
%   seti.dSMeas = 1; % see setMeasPnts.m
%
%   ThetaiOut = minTolIn(ThetaiOut,pdaNv,iOut,FFq,seti)
%
%% Input Arguments
%
% * ThetaiOut   :   inner tolerance of last outer iteration
% * pdaNv       :   numbers of inner iterations stored in a vector
%                   (size 1 x seti.nOut, where nOut is the maximal number of outer iterations)
% * iOut        :   iteration number of outer iterations
% * FFq = FFqMeas   :   scattered field evaluated on receivers positions
%                       (matrix of size seti.measNb x seti.incNb)
%
% Parameters for inner tolerance, see Section "More About".
% We give values of parameters, which worked fine in [2].
%
% * seti.ThetaStart = 0.925
% * seti.ThetaMax = 0.95
% * seti.TolGamma = 0.90
%
% Other parameters in struct seti:
%
% * |seti.tau|    :   Discrepancy principle stops at parameter tau . delta.
% * |seti.delta|  :   relative noise level
% * |seti.FmeasDelta| :   scattered field at receivers positions (i.e. 'the data')
%                         (uScaRX = uTotRX - uIncRX, i.e. total field minus incident field)
%                         (complex matrix of size seti.measNb x seti.incNb)
%
% Some more parameters in struct seti are given in the Example above.
%
% *Optional Input Argument*
%
% * dispDepth   : Depth of displayed messages (0: no, 4 or greater: yes).
%
%% Output Arguments
%
% * ThetaiOut   :   inner tolerance
% 
%% More About
% 
% The primal-dual algorithm is the inner iteration in the computational
% framework. The inner tolerance principle is applied to stop this inner
% iteration within itself.
% 
% This method is described in paragraph 
% "Stopping strategy 2" of Section 4.10 in [2] and
% bases on the inexact stopping rule for a Newton-like method 
% in Chapter 7.5.3 in [1]. We give a brief introduction here.
%
% *Notation*
%
% * $m = \texttt{iOut}$   :   outer iteration number
% * $n$                   :   inner iteration number
% * $\Theta_\mathrm{start} = \texttt{seti.ThetaStart}$
% * $\Theta_{\max} = \texttt{seti.ThetaMax}$
% * $\gamma_\mathrm{tol} = \texttt{seti.TolGamma}$
% * $\Theta_m = \texttt{ThetaiOut}$     :   inner tolerance
%
% *Scheme*
%
% * Compute the quotient of linarized and non-linearized discrepancy, i.e.
%   $\mathrm{dis}_\mathrm{lin}/\mathrm{dis}_\mathrm{nonlin}$.
% * The non-linearized discrepancy is fixed to save computational time.
% * The iteration is stopped 
%   if $\mathrm{dis}_\mathrm{lin}^{(n)}/\mathrm{dis}_\mathrm{nonlin} < \Theta_m$ 
%   for inner tolerance $\Theta_m \in (0,1]$.
%
% *Algorithm: Compute inner tolerance $\Theta_m$*
%
% * Choose $\Theta_\mathrm{start} \in (0,1)$, $\quad$
%  $\Theta_{\max} \in (\Theta_\mathrm{start},1)$, 
%   and $\gamma_\mathrm{tol} \in (0,1]$.
%
% * Auxiliary tolerances for $m = 1,2$:
%
% $\tilde{\Theta}_{m} = \Theta_\mathrm{start}$ 
%
% * Auxiliary tolerances for $m \geq 3$:
%
% With notation: $N_\mathrm{in}^{(m-2)}$ and $N_\mathrm{in}^{(m-1)}$   :   
%  number of steps in inner iteration for the two previous outer
%  iterations.
%
% $\tilde{\Theta}_{m} = 
%   1-\left( 1-\Theta_{m-1} \right) 
%   N_\mathrm{in}^{(m-2)}/N_\mathrm{in}^{(m-1)} 
%   \quad \mathrm{ if\ } N_\mathrm{in}^{(m-1)} \geq N_\mathrm{in}^{(m-2)},$
%
% $\tilde{\Theta}_{m} = 
%   \gamma_\mathrm{tol}\, \Theta_{m-1} \quad \mathrm{ otherwise}.$
%
% * Compute inner tolerance:
%
% $\Theta_m = \Theta_{\max}\, \max \left\{
%   \tau\,\delta/\mathrm{dis}_\mathrm{nonlin}^{(m)},\tilde{\Theta}_{m}
%   \right\}, \quad m \in \bf{N}.$
%
%
%
%% References
% 
% * [1] Andreas Rieder. _Keine Probleme mit Inversen Problemen._ Vieweg, Wiesbaden, 2003.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
% 
%
%% See Also
% * <minPda.html>
% * <pda.html>
%
%
%% Code
%
function ThetaiOut = minTolIn(ThetaiOut,pdaNv,iOut,FFq,seti,varargin)

if nargin == 6
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

seti.useTolInFix = 0; % default: 0
% Do not use option tolInFix because it needs more outer iterations...
if seti.useTolInFix
    ThetaiOut = seti.ThetaStart;
else

    ThetaMax = seti.ThetaMax;
    gamma = seti.TolGamma;

    % Compute help tolerances \hat{Theta}_i = ThetaHati
    if iOut > 2
        if pdaNv(iOut-1) >= pdaNv(iOut-2)
            ThetaHati = 1-pdaNv(iOut-2)/pdaNv(iOut-1)*(1-ThetaiOut);
        else
            ThetaHati = gamma*ThetaiOut;
        end
    else
        ThetaHati = seti.ThetaStart; % start...
    end

    disPri = seti.tau*seti.delta/((normws(FFq-seti.FmeasDelta,seti)/normws(seti.FmeasDelta,seti)));
    ThetaiOut = ThetaMax*max(disPri,ThetaHati); % Rieder p. 265

    %fprintf('        disPri = %1.2g | ThetaHati = %1.2g | ThetaiOut = %1.2g \n',disPri,ThetaHati,ThetaiOut)

    if dispDepth >= 4
        disp('----- inner tolerance principle -----')
        if iOut > 2
            fprintf('r_{n-1} = %g | r_{n-2} = %g',pdaNv(iOut-1),pdaNv(iOut-2));
        end
        fprintf(' | ThetaHati = %g | disPri = %g | ThetaiOut = %g\n',ThetaHati,disPri,ThetaiOut);
        disp('-------------------------------------')
    end
end

end
