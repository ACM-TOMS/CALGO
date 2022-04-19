%% minTolOut
% Outer tolerance principle to stop the primal-dual algorithm.
% 
%% Syntax
%
%   pdaN = minTolOut(iOut,pdaNv,dis,disLin,relDis,seti)
%   pdaN = minTolOut(iOut,pdaNv,dis,disLin,relDis,seti,dispDepth)
%
%% Description
%
% |pdaN = minTolOut(iOut,pdaNv,dis,disLin,relDis,seti)|
% computes the number of inner iterations |pdaN| to stop the primal-dual
% algorithm.
% See "More About" for details.
%
% |pdaN = minTolOut(iOut,pdaNv,dis,disLin,relDis,seti,dispDepth)| does the same but allows to control
% the displayed messages by |dispDepth|.
%
% To use it in context of computational framework set seti.useTolOut = 1.
%
%% Input Arguments
%
% * iOut        :   iteration number of outer iterations
% * pdaNv       :   numbers of inner iterations stored in a vector
%                   (size 1 x seti.nOut, where nOut is the maximal number of outer iterations)
% * dis         :   discrepancy of each outer iteration stored in vector
%                   (size 1 x seti.nOut, where nOut is the maximal number of outer iterations)
% * disLin      :   linearized discrepancy of each outer iteration stored in vector
%                   (size 1 x seti.nOut, where nOut is the maximal number of outer iterations)
% * relDis      :   relative discrepancy =  disLin/dis of each outer iteration stored in vector
%                   (size 1 x seti.nOut, where nOut is the maximal number of outer iterations)
%
% * |seti.pdaNmax|      : Maximal iteration number, if |seti.useTolIn| or
%                         |seti.useTolOut| are 1 (default : 250).
% * |seti.relDisTol|    : Outer tolerance (default: 0.05).
%
% *Optional Input Argument*
%
% * dispDepth   : Depth of displayed messages (0: no, 4 or greater: yes).
%
%% Output Arguments
%
% * pdaN        :   number of inner iterations to stop the primal-dual algorithm
%
%% More About
%
% The primal-dual algorithm is the inner iteration in the computational
% framework. The outer tolerance principle defines a stop index |pdaN|
% for the next application of primal-dual algorithm.
% 
% This method is described in paragraph 
% "Stopping strategy 1" of Section 4.10 in [2].
% We give a brief introduction here.
%
% *Notation*
%
% * $\tau_\mathrm{out} = \texttt{seti.relDisTol}$
%
% *Scheme*
%
% * Choose the number of inner iterations in the $m$-th outer iteration:
%
% $\texttt{pdaN} = N^{(m)}_\mathrm{in}$
%
% * Review this choice: 
% It is "good", if the following quotient is approximately 1:
%
% $\mathrm{dis}_\mathrm{rel}
%  := \mathrm{dis}_\mathrm{lin} / \mathrm{dis}_\mathrm{nonlin}$
%
% (i.e. linearized discrepancy  / non-linearized discrepancy)
%
% If the choice was "good" we increase the number of inner iterations
% (first case), if not we decrease it (second case):
%
% 1st case: $\quad N_\mathrm{in}^{(m+1)} = \lceil \mu_\uparrow\, N_\mathrm{in}^{(m)} \rceil
%  \quad \mathrm{ if\ } \mathrm{dis}_\mathrm{rel} \in
%  (1-\tau_\mathrm{out},1+\tau_\mathrm{out}),$
%
% 2nd case: $\quad N_\mathrm{in}^{(m+1)} =
%  \lfloor \mu_\downarrow\, N_\mathrm{in}^{(m)} \rfloor 
%  \quad \mathrm{ otherwise}.$
%
% Note, that $\mu_\uparrow = 1+\sqrt{1/m}\, \ln(m)$ and 
% $\mu_\downarrow = (\min\{ 1/\mathrm{dis}_\mathrm{rel}, \mathrm{dis}_\mathrm{rel} \})^2$.
%
% * Start with conservative choice: $N_\mathrm{in}^{(1)} = 1$.
% * Limit the number of inner iterations by $\texttt{seti.pdaNmax} = 250$.
%
%
%% References
% 
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <minPda.html>
% * <pda.html>
%
%% Code
%
function pdaN = minTolOut(iOut,pdaNv,dis,disLin,relDis,seti,varargin)

if nargin == 7
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

relDisTol = seti.relDisTol; % tolerance, e.g. 0.01

if iOut > seti.iOutIni+1
    % Adapt the number of inner iterations pdaN.
    
    if dispDepth >= 4
        disp('------------------------------------------')
        fprintf('    pdaN     = %g\n',pdaNv(iOut-1));
        fprintf('dis (nonlin) = %g\n',dis(iOut-1));
        fprintf('dis (lin)    = %g\n',disLin(iOut-1));
        disp('------------------------------------------')
    end

    % Relative difference between linearized and non-linearized
    % discrepancy terms.
    RelDisVal = relDis(iOut-1);
    
    mu1 = 1+sqrt(1./iOut).*log(iOut);
    mu2 = min(1/RelDisVal,RelDisVal)^2;
    pdaNprev = seti.pdaN;
    if dispDepth >= 4
        fprintf('|relDis| = %g\n',RelDisVal);
    end
    if 1-relDisTol < RelDisVal && RelDisVal < 1+relDisTol
        if dispDepth >= 4
            fprintf('%g < |relDis| < %g\n',1-relDisTol,1+relDisTol);
        end
        pdaN = ceil(mu1*pdaNprev); % increase pdaN
    else
        if dispDepth >= 4
            fprintf('|relDis| < %g OR |relDis| > %g\n',1-relDisTol,1+relDisTol);
        end
        pdaN = floor(mu2*pdaNprev); % decrease pdaN
    end
    if pdaN == 0
        pdaN = 1;
    end
    if dispDepth >= 4
        fprintf('pdaN = %g (prev.: %g)\n',pdaN,pdaNprev);
    end
elseif seti.useTolIn == 0
    pdaN = 1; % inner tolerance principle is not used, then do only 1 step.
else
    pdaN = seti.pdaN; % no change
end

if pdaN > seti.pdaNmax
    pdaN = seti.pdaNmax;
end

end
