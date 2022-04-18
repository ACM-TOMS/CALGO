%% setInvType
% Set reconstruction parameters in dependence of seti.invNo 
% (inversion method number).
%
%% Syntax
%
%   seti = setInvType(seti)
%   seti = setInvType(seti,dispDepth)
%
%% Description
% |seti = setInvType(seti)| sets the reconstruction parameters in the
% structure array seti in dependence of the inversion method number
% seti.invNo.
%
% |seti = setInvType(seti,dispDepth)| does the same but allows to control
% the displayed messages by |dispDepth|.
%
% * If no *seti.invNo* was set, set *default 6*.
% * In public version only invNo 6 is available.
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
% * |seti.invNo|  : Number of reconstruction setting (default: 6)
%                   (invNo = 6 is the only one supported in public version)
% * |seti.alpha|  : Regularization parameter for sparsity penalty 
%                   (default: 500).
% * |seti.beta|   : Regularization parameter for total variation penalty
%                   (default: 1E-5).
% * |seti.tau|    : Discrepancy principle stops at parameter tau * delta
%                   (default: 2.5). (Note that delta is the relative noise
%                   level |seti.delta|.)
%
%% Output Arguments
%
% * |seti.inv|        : Type of reconstruction method: 'pda' or 'shrinkage'
% * |seti.useWavelet| : Using wavelets (1) or not (0)
%                       (Wavelets are not available in public version)
% * |seti.tF|     : Terms belonging to functional F
%                   (e.g. |seti.tF = {'fd','fg'}|; See table below.)
% * |seti.tG|     : Terms belonging to functional G
%                   (e.g. |seti.tG = {'fs','fp'}|; See table below.)
% * |seti.pNorm|  : Index p for p-Schatten-Norm of measurements to define
%                   discrepancy.
% * _More about pNorm:_   default: pNorm = 2, i.e. essentially Frobenius-Norm
% * _More about pNorm:_   in case of inv = 'pda' with invNo 3 to 6
%                         set pNorm = 2 (automatically) because it is required
% * |seti.qNorm|      : Index q for q-Norm of contrasts to define sparsity penalty
% * _More about qNorm_: Standard: qNorm = 1 -> l^1 minimization.
% * _More about qNorm_: in case of inv = 'pda' set qNorm = 1 (automatically)
%                       because it is required.
% * |seti.p|      : Exponent in discrepancy-term $1/p ||...||^p_2$ with $p = 1$ or $2$
%                   (used e.g. in <setFuncsPda.html>)
%                   (default: |seti.p = 2|).
%
% Note: *p ~= pNorm* in general!
%                   (but: in case of inv = 'shrinkage':
%                   pNorm and qNorm are used as exponents too.)
%
%% More About
%
% *Shrinkage*
%
% * You can influence method changing
% seti.pNorm (default 2), seti.qNorm (default 1).
% * Parameter seti.p has no influence in shrinkage.
%   (In case of shrinkage(!) exponents are always seti.pNorm and seti.qNorm.)
%
% *pda (primal-dual algorithm)*
%
% * In pda must be seti.pNorm = 2, seti.qNorm = 1 (except for invNo = 7).
%   (In cases of invNo = 3:1:6 this is checked below.)
% * In case of invNo = 7 seti.pNorm must be the same as seti.p. (This case
%   is not restricted to seti.qNorm = 1).
% * Parameter seti.p can be changed only in pda with wavelets.
%   So, you have invNo 4 and 5 (without wavelets it's fixed p = 2).
% * Details of used terms (fd, fs, ...) and their belonging to term F or G
%   are inside <pda.html>.
%
% *invNo: inversion number: 1--8*
%
% * d: discrepancy
% * s: sparsity
% * g: total variation (g for gradient...)
% * p: physical bounds
%
%
%  No. d    s    g   p  (explanation)
%  1) fd   fs   --- fp  (shrinkage without wavelets)
%  2) |fd  fsw  --- |fp (shrinkage with wavelets)
%  3) fd   fs   --- fp  (pda without wavelets)
%  4) fdw1 fsw  --- |fp (pda with wavelets, exponent p = 1)
%  5) fdw2 |fsw --- |fp (pda with wavelets, exponent p = 2)
%  6) fd   fs   fg  fp  (pda without wavelets)
%  7) fdw3 fsw  --- |fp (pda with wavelets, exponent p = p)
%  8) fd   fs   fg  fp2 (pda without wavelets; like 6 but with individual 
%                        physical bounds for background)
%
%
% *Minimization functional in case of invNo = 6 (default case)*
%
% $$ \min_{h \in X} 
%      \underbrace{
%          \frac{1}{2}\|\mathcal{F}'(q)[h]+\mathcal{F}(q) - F_\mathrm{meas}^\delta\|_\mathrm{F}^2
%          }_{=: f_\mathrm{dis}(h),\ \mathrm{discrepancy\ (linearized\ problem)}}
%    + \underbrace{\alpha \|q+h\|_\mathrm{1}}_{=: f_\mathrm{spa}(h),\ \mathrm{sparsity\ penalty}}
%    + \underbrace{\beta \| \nabla (q+h) \|_\mathrm{1}}_{=: f_\mathrm{tv}(h),\ \mathrm{total\ variation\ penalty}}
%    + \underbrace{
%          \delta_{[a,b]}( \mathrm{Re}(q+h) ) +
%          \delta_{[c,d]}(\mathrm{Im}(q+h) )}_{=: f_\mathrm{phy}(h),\ \mathrm{penalty\ for\ physical\ bounds}
%          }.
%    $$
%
% *Minimization functional in a more general case*
%
% Notation:
%
% * $p = \texttt{seti.p}$
% * $P = \texttt{seti.pNorm}$
% * $Q = \texttt{seti.qNorm}$
%
% $$ \min_{h \in X} 
%      \underbrace{
%          \frac{1}{p}\|\mathcal{F}'(q)[h]+\mathcal{F}(q) -
%          F_\mathrm{meas}^\delta\|_\mathrm{Schatten\ P-Norm}^p
%          }_{\mathrm{discrepancy\ (linearized\ problem)}}
%    + \underbrace{\alpha \|q+h\|_Q}_{\mathrm{sparsity\ penalty}}
%    + \underbrace{\beta \| \nabla (q+h) \|_\mathrm{1}}_{\mathrm{total\ variation\ penalty}}
%    + \underbrace{
%          \delta_{[a,b]}( \mathrm{Re}(q+h) ) +
%          \delta_{[c,d]}(\mathrm{Im}(q+h) )}_{\mathrm{penalty\ for\ physical\ bounds}
%          }.
%    $$
%
%% See Also
%
% * <setRecon.html>
% * <pda.html>
%
%% Code
%
function seti = setInvType(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

seti = checkfield(seti,'invNo',6,dispDepth); 

% idea for other name "invset" (inv setting...) (or name in suffix of dir...)

% "working parameters", i.e. some working parameters in case of contrast twoCornersOneBall2D

switch seti.invNo
    case 1 % shrinkage, no wavelets
        %shrinkWavNo
        seti.inv = 'shrinkage';
        seti.useWavelet = 0;
        
        % working parameters
        seti = checkfield(seti,'alpha',0.005,dispDepth); % was 0.002 until 201704
        
    case 2 % shrinkage, with wavelets
        %shrinkWavYes
        seti.inv = 'shrinkage';
        seti.useWavelet = 1;
        
        % working parameters
        seti = checkfield(seti,'alpha',0.02,dispDepth);
       
    case 3 % pda, without wavelets, without fg
        %pdaWavNoGradNo
        seti.inv = 'pda';
        seti.useWavelet = 0;
        % without wavelets, without fg term
        % (and p = 2, but p is not set, because ^2 in Code)
        % F(x) = fd = 1/2 || x + FF(q) - FmeasDelta||WS,2^2
        % K(x) = FF'(q)[x]
        % So: F(K(x)) = 1/2 || FF'(q)[x] + FF(q) - FmeasDelta||WS,2^2
        % G(x) = fs + fp = \alpha ||Re(q+x)||_1 + \alpha ||Im(q+x)|| + delta_[0,c](Im(q+x) + delta_[a,b](Re(q+x))
        seti.tF = {'fd'};
        seti.tG = {'fs','fp'}; % cell array of strings

        % working parameters
        seti = checkfield(seti,'alpha',2,dispDepth);
        
    case 4 % pda, with wavelets and p = 1
        %pdaWavYesp1
        seti.inv = 'pda';
        seti.useWavelet = 1;
        seti.p = 1; % Attention: p ~= pNorm in general!
        seti.tF = {'fdw1','fsw'};
        seti.tG = {'fp'};
        
        % working parameters
        seti = checkfield(seti,'alpha',0.05,dispDepth); 

    case 5 % pda, with wavelets and p = 2
        %pdaWavYesp2
        seti.inv = 'pda';
        seti.useWavelet = 1;
        seti.p = 2; % Attention: p ~= pNorm in general!
        seti.tF = {'fdw2','fsw'};
        seti.tG = {'fp'};
        
        % working parameters
        seti = checkfield(seti,'alpha',0.3,dispDepth); 

    case 6 % pda, without wavelets, with fg term
        % (and p = 2, but p is not set, because ^2 in Code)
        %pdaWavNoGradYes
        seti.inv = 'pda';
        seti.useWavelet = 0;
        seti.tF = {'fd','fg'};
        seti.tG = {'fs','fp'};
        
        % working parameters
        seti = checkfield(seti,'alpha',500,dispDepth);
        seti = checkfield(seti,'beta',1E-5,dispDepth);
        seti = checkfield(seti,'tau',2.5,dispDepth);
        
    case 7 % pda with wavelets: p-Norm in discrepancy, q-Norm in penalty
        % (This case requires p = P, i.e. seti.p = seti.pNorm.)
        %
        % -- Overview of the formulas:
        %
        % F(x) = fd + fs with
        % fd = 1/P || x + FF(qW) - FmeasDelta||WS,P^P
        % fs = (\alpha/Q) ||Re(q+x)||_Q^Q + (\alpha/Q) ||Im(q+x)||_Q^Q
        % (To keep the notation simple waveletes are omitted in the notation of fs.)
        % K(x) = [Kdis(x), Kspa(x)] with Kdis(x) = FF'(qW)[x]
        % So: F(K(x)) = 1/P || FF'(qW)[x] + FF(qW) - FmeasDelta||WS,P^P
        %               + sparsity term
        % G(x) = fp = delta_[0,c](Im(q+x) + delta_[a,b](Re(q+x))

        seti.inv = 'pda';
        seti.useWavelet = 1;
        seti.tF = {'fdw3','fsw'};
        seti.tG = {'fp'};
        
        % working parameters
        seti = checkfield(seti,'alpha',0.03,dispDepth);
        
    case 8 % (pda without wavelets; like 6 but with individual physical bounds for background)
        seti.inv = 'pda';
        seti.useWavelet = 0;
        seti.tF = {'fd','fg'};
        seti.tG = {'fs','fp2'};
        
        % working parameters
        seti = checkfield(seti,'alpha',500,dispDepth);
        seti = checkfield(seti,'beta',1E-5,dispDepth);
        seti = checkfield(seti,'tau',2.5,dispDepth);
        
end

%%
% *Check pNorm and qNorm compatibility to cases*
%
if ((3 <= seti.invNo && seti.invNo <= 6) || seti.invNo == 8) && (~isfield(seti,'pNorm') || seti.pNorm ~= 2)
    seti.pNorm = 2;
    if dispDepth >= 1
        disp('   PDA cases with invNo 3 to 6 and 8 requires pNorm = 2. Set it.');
    end
end

if strcmp(seti.inv,'pda') && (~isfield(seti,'qNorm') || seti.qNorm ~= 1)
    seti.qNorm = 1;
    setmessage(seti,'qNorm',dispDepth);
    if dispDepth >= 1
        disp('   (To use pda, you have to use parameter qNorm = 1.)');
    end
end


end
