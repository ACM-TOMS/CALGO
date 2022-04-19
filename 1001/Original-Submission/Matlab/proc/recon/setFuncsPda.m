%% setFuncsPda
%
% Set functions required in primal-dual algorithm (pda.m).
%
%% Syntax
%
%   seti = setFuncsPda(seti)
%
%% Description
% |seti = setFuncsPda(seti)| sets functions in structure array |seti|,
% that are required in primal-dual algorithm, see <pda.html>. 
% Therefore this function is called in setRecon.m, see <setRecon.html>.
%
%% Input Arguments
%
% * seti    :   structure array
%
%% Output Arguments
%
% * seti    :   structure array
%
% For the fields of the structrual array |seti| see the section "More About".
%
%% More About
%
% We only give a overview of the functions. See also <pda.html>.
% Detailed information are in Section 4.5 in [1].
%
% Note that grid scaling is respected in the code, but for 
% simplicitiy is not presented in the formulas.
%
% *Auxiliary quantities*
%
% For $\texttt{seti.S} = T_{\bf{C}\to\bf{R}^2}$ and $\texttt{seti.T} = T_{\bf{R}^2\to\bf{C}}$, see "More About" in <pda.html>.
%
% * $\texttt{seti.Kd} = K_\mathrm{dis} :=
%   T_{\bf{C}\to\bf{R}^2} [\mathcal{F}'(T_{\bf{R}^2\to\bf{C}}\ q)]T_{\bf{R}^2\to\bf{C}}$,
% * $\texttt{seti.vd} = v_\mathrm{dis} := 
%   T_{\bf{C}\to\bf{R}^2}(\mathcal{F}(T_{\bf{R}^2\to\bf{C}}\ q)-F_\mathrm{meas}^\delta)$,
% * $\texttt{seti.Kg} = K_\mathrm{tv}  := 
%   \beta\,T_{\bf{C}\to\bf{R}^2}\nabla T_{\bf{R}^2\to\bf{C}}$
% * $\texttt{seti.vg} = v_\mathrm{tv} :=
%   \beta\, T_{\bf{C}\to\bf{R}^2} \nabla T_{\bf{R}^2\to\bf{C}}\ q$.
%
% Note that the names |Kg| and |vg| were choosen because the appearance of the gradient.
%
% *Parts of functional to minimize*
%
% * $\texttt{seti.fd} = f_\mathrm{dis}(h) := 
%   \frac{1}{2} \|K_\mathrm{dis} h+  v_\mathrm{dis} \|_{\mathrm{dis},\bf{R}}^2$
%   $\quad$ is the discrepancy (linearized problem),
% * $\texttt{seti.fs} = f_\mathrm{spa}(h) := 
%   \alpha \| h + q\|_{\mathrm{spa}, \bf{R}}$
%   $\quad$ is the sparsity penalty,
% * $\texttt{seti.fg} = f_\mathrm{tv}(h) :=
%   \|K_\mathrm{tv} h + v_\mathrm{tv} \|_{\mathrm{tv},\bf{R}}$
%   $\quad$ is the TV-penalty,
% * $\texttt{seti.fp} = f_\mathrm{phy}(h) :=
%   \delta_{[a,b,c,d]}(h + q)$
%   $\quad$ is the penalty for physical bounds.
%
% Note that $\delta_{[a,b,c,d]}(x)$ is a indicator function, i.e. is 0, 
% if all entries of the vector $\mathrm{real}(x)$ are between $a$ and $b$
% and all entries of the vector $\mathrm{imag}(x)$ are between $c$ and $d$;
% and is $\infty$, otherwise.
%
% *Definition of parts of Tikhonov functional*
%
% The Tikhonov functional is |MT = M1 + M2|. 
% Function $F$ is $M_1$, function $G$ is $M_2$ (with $F$ and $G$ from PDA)
% This is done in function |minPda| (<minPda.html>) when function 
% |pda| (<pda.html>) is called.
%
% *Further comments (in particular, interesting in not public version)*
%
% * In case of inversion by pda, it is always (automatically) set 
%   |seti.pNorm = 2| and |seti.qNorm = 1| (except for invNo = 7), see
%   <setInvType.html>.
% * Exponent p of the discrepancy term, see <setInvType.html>, is usually
%   |seti.p = 2|.
%   In case of shrinkage or pda with wavelets it may differ, see
%   <setInvType.html>.
% * |vs| and |vp| in term $G$: If you choose |vs| and |vp| manually and 
%   |fs| and |fp| belong to term $G$, then |vs = vp| is required.
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <setRecon.html>
% * <pda.html>
% * <setInvType.html>
%
%% Code: setFuncsPda
%
function seti = setFuncsPda(seti)

% Definition of parts of function to minimize (f = minf)
% parts for reconstruction process
seti = setminfParts(seti);
if seti.useWavelet == 1 % not available in public version
    seti = setminfPartsWavelet(seti);
end

end

%% Code: subfunction: potPotNum
% Numerical potential pot function.
%
function [d1,d2] = potPotNum(seti)
infty = 1E5;
dp = @(a,b,x) 0.* and((a<=x),(x<=b))+infty.*or((x<a),(b<x)); % potential pot function
d =  @(a,b,x) max(dp(a,b,x)); % d = 0 if ALL entries are between a and b, otherwise d = \infty
% max because d = \infty if ONE entry of x \notin [a,b]
clear infty

l = seti.physBounds;
reMin = l(1); reMax = l(2); imMin = l(3); imMax = l(4);
tol = seti.tolDelta;
d1 = @(x) d(reMin-tol,reMax+tol,x); % limits for real part
d2 = @(x) d(imMin-tol,imMax+tol,x); % limits for imag part
% tol important because e.g. (min(real(qROI))+1) = -2.2204e-16
clear l a b c
end

%% Code: subfunction: setminfParts
%
% * $K: X \to Y$
% * Y: depends on function |Ki| (|Ki| i.e. |Kd|, |Kg|, ...)
% * X: xnRVD (this means real vector downscaled):
%    size of elements: 2*seti.nInv^seti.dim
% * |hs| was |xnRVU|
% * |q| is standard in qCVU (is fixed in pda), 
%   but you can use qCVD etc. too, if you have them as argument

function seti = setminfParts(seti)

%%
% *f_d*
%
% * Remember that |DFFq = @(xnRVD) JA*diag(seti.GU(seti.T(xnRVD)))*JB;|
% * The difference to DFFq is that finally |seti.S| is used.
%
if strcmp(seti.model,'helmholtz2D') || strcmp(seti.model,'helmholtz3D')
    seti.Kd = @(xnRVD,JA,JB) seti.S(JA*diagsparse(seti.GU(seti.T(xnRVD)))*JB); % K_d(h) = S(\F'(q)[GU(T(h_S))]) %real
elseif strcmp(seti.model,'helmholtzHMode2D')
    seti.Kd = @(xnRVD,JA,JB) seti.S( JA(:,:,1)*diagsparse(seti.GU(seti.T(xnRVD)))*JB(:,:,1) ...
        + JA(:,:,2)*diagsparse(seti.GU(seti.T(xnRVD)))*JB(:,:,2) ); % K_d(h) = S(\F'(q)[GU(T(h_S))]) %real
else
    fprintf('Error in setFuncs - pda not implemented for model %s.\n',  seti.model)
end
% Finally: Kd(x): R^{2 nInv^dim} -> R^{2 measNb x incNb}
% Kd = \F'(q) (from mimo \F'(q): JA*diag(h)*JB = \F'(q)[h])
seti.KdAdj = @(yd,JA,JB) seti.S(seti.GD(ADFFqFast(seti.T(yd),JA,JB,seti)));
seti.vd = @(FFq) seti.S(FFq-seti.FmeasDelta); % 2* measNb x incNb
seti.fd = @(DFFqh,FFq) (1/2)*normws(seti.S(DFFqh)+seti.vd(FFq),seti)^2;

%%
% *f_s*
%
seti.vs = @(qCVU) seti.S(seti.GD(qCVU));
seti.fs = @(qCVU,xnRVD) seti.alpha*normroi(xnRVD+seti.vs(qCVU),seti);

%%
% *f_g*
%
seti.Kg = @(xnRVD) seti.beta*gradientNeumann(xnRVD,seti.hInv,seti.nInv,seti.GInv,seti); % K_g
%seti.KgAdj = @(yg) -seti.beta*div(yg,seti); % Note that this is wrong(!) (although analytically correct)
seti.KgAdj = @(yg) seti.beta*gradientNeumannAdj(yg,seti); %yg is 2*dim x nInv x nInv (x nInv)
seti.vg = @(qCVU) seti.Kg(seti.S(seti.GD(qCVU)));
seti.fg = @(qCVU,xnRVD) normTVinv1(seti.Kg(xnRVD)+seti.vg(qCVU),seti); %f_g(x) = ||beta*\nabla x + v_g)||_\tv (gradient sparsity)

%%
% *f_p*
%
% * Note that seti.vp = seti.vs is important because of prox Gsp
%   (prox Gsp does mean that the term G consists of fs and fp.
%   To compute the proximal mapping of G the equality of vp and vs is
%   required!)
%
[d1,d2] = potPotNum(seti);
%seti.fp = @(q,hz) d1(real(q)+real(hz)) + d2(imag(q)+imag(hz));
% use seti.vs in fp:
%seti.fp = @(q,hz) d1(seti.R(seti.vs(q))+seti.R(seti.vs(hz))) + d2(seti.I(seti.vs(q))+seti.I(seti.vs(hz)));
seti.fp = @(qCVU,xnRVD) d1(real(seti.GD(qCVU))+seti.R(xnRVD)) + d2(imag(seti.GD(qCVU))+seti.I(xnRVD));

seti.fp2 = seti.fp; % In fact it is not the same function, but because the pixels in the background are zero the result is the same.

end

%% Code: subfunction: setminfPartsWavelet
% Definitions of functions in case of wavelets.
% (This is not supported in the public version.)
%
function seti = setminfPartsWavelet(seti)

% f_dw
% seti.Kdw = seti.Kd;
%seti.KdwAdj = seti.KdAdj;
seti.vdw = @(FFq) seti.S(FFq-seti.FmeasDelta);
%seti.fdw = @(DFFqh,FFq) 1/seti.p*normws(seti.S(DFFqh)+seti.vdw(FFq),seti)^seti.p;

seti.fdw1 = @(DFFqh,FFq) normws(seti.S(DFFqh)+seti.vdw(FFq),seti); % seti.p = 1
seti.fdw2 = @(DFFqh,FFq) 1/2*normws(seti.S(DFFqh)+seti.vdw(FFq),seti)^2; % seti.p = 2
seti.fdw3 = @(DFFqh,FFq) 1/seti.pNorm*normLp(DFFqh+FFq,seti)^seti.pNorm; % seti.qNorm ~= 2
%fswPart = @(q,hz) seti.wW(q+hz);
%seti.fsw = @(q,hz) normroi(seti.alpha*[real(fswPart(q,hz)); imag(fswPart(q,hz))],seti);

% f_sw
if isfield(seti,'wavIsom') && strcmp(seti.wavIsom,'W1p') && isfield(seti,'omegaW1p')
    seti.Ksw = @(xnRVD) seti.alpha*seti.S(seti.wWi(seti.omegaW1p.*seti.T(xnRVD)));
    seti.KswAdj = @(ys) seti.alpha*seti.S(seti.omegaW1p.*seti.wWstari(seti.T(ys)));
    %seti.KswAdj = @(ys) seti.alpha*[seti.wWstar(seti.R(ys)); seti.wWstar(seti.I(ys))];
    %seti.vsw = @(q) seti.alpha*[real(wWq(q)); imag(wWq(q))];
    seti.vsw = @(qCVU) seti.alpha*seti.S(seti.wWi(seti.omegaW1p.*seti.GD(qCVU))); % down scaled
    seti.fsw = @(qCVU,xnRVD) norminv(seti.vsw(seti.GU(seti.T(xnRVD)))+seti.vsw(qCVU),seti);
    % T(xnRVD) is upscaled and in seti.vsw down scaled...
    % using qCVD (down scaled) instead of qCVU would be useful...
    % In theory normroi(seti.S(...),...), but it does not matter
else
    seti.Ksw = @(xnRVD) seti.alpha*seti.S(seti.wWi(seti.T(xnRVD)));
    seti.KswAdj = @(ys) seti.alpha*seti.S(seti.wWstari(seti.T(ys)));
    seti.vsw = @(qCVU) seti.alpha*seti.S(seti.wWi(seti.GD(qCVU))); % down scaled
    seti.fsw = @(qCVU,xnRVD) norminv(seti.vsw(seti.GU(seti.T(xnRVD)))+seti.vsw(qCVU),seti);
end


end
