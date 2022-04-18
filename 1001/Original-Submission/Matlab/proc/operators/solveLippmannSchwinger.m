%% solveLippmannSchwinger
% Solve Lippmann-Schwinger equation.
%
%% Syntax
%
%   v = solveLippmannSchwinger(Vq, f, seti, method, Vqf, VqfM)
%
%% Description
% |v = solveLippmannSchwinger(Vq, f, seti, method, Vqf, VqfM)|
% solves the Lippmann-Schwinger equation: computes v such that v - Vq(v) = f
% for a function Vq.
%
%% Input Arguments
%
% * Vq       : function in Lippmann-Schwinger equation |v - Vq(v) = f|
% * f        : f in Lippmann-Schwinger equation |v - Vq(v) = f|
% * seti.tol : tolerance to stop the solving process
%
% _Optional input arguments:_
%
% * method  : solving methods: 'GMRES', 'preconditionedGMRES', 
%             'GMRESkelley' (default), 
%             'twoGrid' (not supported in public version)
% * VqM, fM : equivalents of Vq, f on the coarse grid, 
%             needed for method='twoGrid', which is
%             not supported in the public version.
%             Note that twoGrid is only used if seti.mCD > 0.
%
%% Output Arguments
%
% * v       : solution of the Lippmann-Schwinger equation:
%             |v| solves |v - Vq(v) = f|.
%             (In context of this framework $v$ is the scattered field on
%             ROI, i.e. |uScattROI|)
%
%% More About
%
% * The Lippmann-Schwinger equation is an integral equation.
% * For the Lippmann-Schwinger equation see [1] or [2, Sec. 3.1].
%
% *Lippmann-Schwinger equation in the code*
%
% $v$ is a solution of the Lippmann-Schwinger equation:
%
% $v - \texttt{Vq}(v) = f$.
%
% *Lippmann-Schwinger equation in the context of scattering*
%
% In the context of scattering we use the Lippmann-Schwinger equation
% to compute the scattered field $u^s = \texttt{uScattROI}$ on ROI, see [1] or [2, Sec. 3.1]:
%
% $u^s$ solves the Lippman-Schwinger equation:
%
% $u^s - V(q.*u^s) = V(q.*u^i)$
%
% with incident field $u^i$ on ROI, contrast $q$ on ROI and volume
% potential $V$ defined by
%
% $(Vf)(x):= k^2 \int_\mathrm{ROI} \Phi(x-y) f(y)\,\mathrm{d}y, \quad x\in \bf{R}^d$
%
% with dimension $d$ and fundamental solution $\Phi$ 
% (for a definition see <setIncField.html>, Section "More About").
%
% Transferred to code this results in, see <simo.html>:
%
%   uScattROI = solveLippmannSchwinger(@(x) V(QU(x)), V(QU(uIncROI)), seti); 
%
% with
%
% * $\texttt{Vq} = \texttt{@(x) V(QU(x))} = \texttt{@(x)} V\texttt{(qROI.*x)}$,
% * $\texttt{f} = \texttt{V(QU(uIncROI))} = V(q.*u^i)$.
% * Note that $V$ is the volume potential defined in <intOpsFuncs.html> by
%   |V  = @(x) seti.k^2.*helmholtz2Dr2r(x, seti);|
%
% *Notation*
%
% * $\texttt{uScattROI} = u^s$ (scattered field on ROI)
% * $\texttt{uIncROI} = u^i$ (incident field on ROI)
%
% 
%% References
%
% * [1] David Colton and Rainer Kress. _Inverse Acoustic and Electromagnetic Scattering Theory_. Springer, New York, 2013.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <simo.html>
% * <intOpsFuncs.html>
% * <setIncField.html>
%
%% Code
% *Code: solveLippmannSchwinger*

function v = solveLippmannSchwinger(Vq, f, seti, method, Vqf, VqfM)

tol = seti.tol;

if ~exist('method','var') || (strcmpi(method,'twoGrid') && seti.mCD <= 0) 
    method = 'GMRESkelley';
end

A = @(x) x - Vq(x);
if strcmpi(method, 'GMRES')
    [v,flag] = gmres(A , f(:), 15, tol, 300); %#ok<ASGLU>
elseif strcmpi(method, 'preconditionedGMRES')
    D = @(x) precondition(Vq, x, griddata.dastruCoar);
    [v, flag] = gmres(A, f(:), 15, tol, 300, D); %#ok<ASGLU>
elseif strcmpi(method, 'GMRESkelley')
    v = gmresKelley(f(:), f(:), A , [tol, 100, 1]);
elseif strcmpi(method, 'twoGrid')
    v = solveLippmannSchwingerTwoGrid(seti,Vqf,VqfM);
end

end

%%
% *Code: subfunction: precondition*

function x = precondition(Vq,x,griddata)

% preconditioner for gmres
Nfine = griddata.NFine;
N     = griddata.N;

x = downsample(x,Nfine,N);
A = @(x) x - Vq(x,griddata);
% [x, flag] = gmres(A ,x(:), 10, 1E-6, 300); %#ok<NASGU>
x = gmresKelley(x(:), x(:), A , [tol ,100,1]);
x = reshape(x, N, N);
x = interpolate(x, N, Nfine);
end
