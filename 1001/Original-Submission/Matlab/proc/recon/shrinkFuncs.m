%% shrinkFuncs
% Defines extended soft-shrinkage functions (for real and imaginary part).
%
%% Syntax
%
%  seti = shrinkFuncs(seti)
%
%% Description
% |seti = shrinkFuncs(seti)| defines extended soft-shrinkage functions for
% real and imaginary part in |seti.shrkRe(x,alpha)| and
% |seti.shrkIm(x,alpha)|.
% 
% * In public version only the case |seti.qNorm = 1| is supported.
%
%% Example
%
% *Example 1 (without physical bounds, i.e. soft-shrinkage)*
%
%   seti.qNorm = 1;
%   seti.physBounds = [-inf, inf, -inf, +inf];
%   seti = shrinkFuncs(seti);
%   alpha = 0.5;
%   t = -4:0.1:4;
%   figure(101);
%   plot(t,seti.shrkRe(t,alpha));
%
% *Example 2 (with physical bounds, i.e. extended soft-shrinkage)*
%
%   seti.qNorm = 1;
%   seti.physBounds = [-1, 3, 0, 3];
%   seti = shrinkFuncs(seti);
%   alpha = 0.5;
%   t = -2:0.1:4;
%   figure(102);
%   plot(t,seti.shrkRe(t,alpha));
%
%
%% Input Arguments
%
% * seti    :   structure array
%
% * seti.qNorm  :   in sparsity term $f_\mathrm{spa}$ a norm with index 
%                   |seti.qNorm = 1| is used 
%                   (other values are not supported in public version).
%                   (Index q = 1 does essentially (up to some factor) mean $l^1$ minimization)
% * seti.physBounds     :   contains the physical bounds (min and max)
%                           for real and imaginary part of the contrast as
%                           4-vector with structure 
%                           |[reMin reMax imMin imMax]|.
%
%% Output Arguments
%
% * seti.shrkRe = @(x,alpha)    :   extended soft-shrinkage function for real part
% * seti.shrkIm = @(x,alpha)    :   extended soft-shrinkage function for imaginary part
%
%% More About
%
% *Soft-shrinkage and extended soft-shrinkage operator*
%
% <<../extGraph/shrinkFuncs_proxg_05.png>>
%
% * red continuous  :   soft-shrinkage
% * blue broken     :   extended soft-shrinkage
% 
% The *soft-shrinkage operator*, see [1], is defined for real-valued $x$ and $a \geq 0$ by
%
% * $\mathcal{S}(x, \kappa) := x + \kappa \quad \mathrm{ if\ } x \leq -\kappa$,
% * $\mathcal{S}(x, \kappa) := 0 \quad \mathrm{ if\ } x \in (-\kappa,+\kappa)$,
% * $\mathcal{S}(x, \kappa) := x - \kappa \quad \mathrm{ if\ } x \geq \kappa$.
%
% Note that $\texttt{alpha}$ instead of $\kappa$ is used in the code below.
%
% The *interval projection operator*, see Section 4.7 in [2], 
% is defined for real-valued $x$ and $a\geq b$ by
%
% * $\mathcal{I}_{[a,b]}(x) :=  a \quad \mathrm{ if\ } x < a$, 
% * $\mathcal{I}_{[a,b]}(x) :=  x \quad \mathrm{ if\ } x \in [a,b]$, 
% * $\mathcal{I}_{[a,b]}(x) :=  b \quad \mathrm{ if\ } x > b$.
%
% The *extended soft-shrinkage*, see Section 4.7 in [2], defined as
%
% $\mathcal{I}_{[a,b]}(\mathcal{S}(x,\kappa))$.
%
% *Some further notes*
%
% * Because we take into account physical bounds we call it 
% _extended_ soft-shrinkage instead of the well-known soft-shrinkage.
% * |shrinkFuncs| was originally written for reconstruction process by 
%   shrinkage, but the defined |seti.shrkRe| and |seti.shrkIm| are also 
%   needed in <pda.m>.
% * seti.physBounds : When shrinkage with wavelets is used, no bounds will
%                     be set (physBounds = [-inf,inf,-inf,inf]).
%
%% References
% * [1] Ingrid Daubechies, Michel Defrise, and Christine De Mol. 
%   An iterative thresholding algorithm for linear inverse problems with a sparsity constraint. 
%   _Communications on Pure and Applied Mathematics_, 57(11):1413-1457, 2004.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% Code
%
function seti = shrinkFuncs(seti)
% defines extended soft-shrinkage functions (for real and imaginary part)

reMin = seti.physBounds(1);
reMax = seti.physBounds(2);
imMin = seti.physBounds(3);
imMax = seti.physBounds(4);

if seti.qNorm == 1
    % shrinkage func for real part
    seti.shrkRe = @(x,alpha) max(min(max((abs(x)-alpha),0).*sign(x), reMax), reMin);
    % shrinkage func for imag part
    seti.shrkIm = @(x,alpha) max(min(max((abs(x)-alpha),0).*sign(x), imMax), imMin);
else % shrinkFuncComp is not available in public code
    seti.shrkRe = @(x,alpha) max(min( shrinkFuncComp(x,alpha,seti.qNorm, 'newton'), reMax), reMin);
    seti.shrkIm = @(x,alpha) max(min( shrinkFuncComp(x,alpha,seti.qNorm, 'newton'), imMax), imMin);
end

end
