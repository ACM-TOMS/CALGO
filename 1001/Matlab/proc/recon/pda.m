%% pda
% primal-dual algorithm in context of the 
% Computational Framework for Inverse Medium Problem in Scattering.
%
%% Syntax
%   [hsolCVU,pdaStopInd,FGval,Fval,Gval,relLinDis,disLinInPda,errInPda,minf] = pda(iOut,qCVU,ThetaiOut,seti,dispDepth)
%
%% Description
% |[hsolCVU,pdaStopInd,FGval,Fval,Gval,relLinDis,disLinInPda,errInPda,minf] 
% = pda(iOut,qCVU,ThetaiOut,seti)|
% computes the update |hsolCVU| of the linearized minimization problem.
%
%
%% Input Arguments
%
% Most important input arguments (there are more...)
%
% * iOut        :   number of (current) outer iteration
% * qCVU        :   reconstructed contrast q as vector of size
%                   seti.nROI^seti.dim x 1.
%                   (CVU does mean complex, vector, upscaled)
% * ThetaiOut   :   inner tolerance, see <minTolIn.html>
%                   (only required if inner tolerance principle is used by
%                   seti.useTolIn = 1; otherwise it is set to 0.)
% * seti        :   struct seti
% * dispDepth   :   depth of displayed messages
%                   (greater or equal 3 to see aything in this file).
%                   (3 less, 4 default, 5 additional information).
%
% *Some of the fields in struct seti*
%
% Note that several fields in struct seti are necessary to run pda.
%
% The routine |pda| is an internal one that needs the specific environment created in
% the package.
%
% * seti.dim    :   dimension of the problem (2 or 3)
% * seti.nROI   :   discretization points for each dimension
%                   of region of interest (ROI) (in samples)
%
% * seti.pdaN   :   number of inner iterations (PDA) (is called nPda)
% * seti.pdaStepsize    :   method to choose primal and dual stepsizes
%                           (|'fix'| (default) or |'adaptive'|), 
%                           see <pdaChoosingStepsizes.html>.
% * seti.vartheta       :   parameter in (0,1) used in case of adaptive stepsizes,
%                           see <pdaChoosingStepsizes.html>.
%
%
%% Output Arguments
%
% * hsolCVU     :   solution of the update h
%                   (at the end of pda: hsolCVU = seti.GU(seti.T(xnRVD)))
% * pdaStopInd  :   iteration, where primal-dual agorithm was stopped
%                   (iPda = 1:nPda, but maybe is stopped earlier becaue
%                   usage of seti.useTolIn or seti.useTolOut)
%                   (see <minTolIn.html>, <minTolOut.html>, <minPda.html>).
% * FGval       :   stores values of $F(Kh) + G(h)$ (vector of size nPda x 1)
% * Fval        :   stores values of $F(Kh)$ (vector of size nPda x 1)
% * Gval        :   stores values of $G(h)$ (vector of size nPda x 1)
% * relLinDis   :   quotient: relLinDis = disLin / dis 
%                   (vector of size nPda x 1)
% * disLinInPda :   discrepancy of linearized problem for each inner iteration step
%                   (vector of size nPda x 1)
% * errInPda    :   relative error of the reconstructed contrast qROI
%                   (vector of size nPda x 1)
% * minf        :   struct with parts of the minimization functional
%                   (Tikhonov functional)
%
% * minf.fd     :   discrepancy of linearized problem
%                   (vector of size nPda x 1)
% * minf.fs     :   sparsity penalty
%                   (vector of size nPda x 1)
% * minf.fg     :   total variation penalty
%                   (vector of size nPda x 1)
% * minf.fp     :   penalty for physical bounds
%                   (vector of size nPda x 1)
%
% See "More About" for formulas.
%
%
%% More About
%
% *Minimization functional*
%
% The function pda minimizes the *Tikhonov functional* of the linearized
% problem, i.e. 
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
% *Next step*
%
% The update $q := q + h$ is done in <minPda.html>, see also "More About" 
% in <start.html> to see the connection.
%
% *Splitting the minimization problem*
%
% In pda we formulate the minimization problem as:
% 
% $$ \min_{h \in X} F(Kh) + G(h)$$
% 
% with
% 
% $F(Kh) = f_\mathrm{dis}(h) + f_\mathrm{tv}(h)$,
% 
% $G(h)  = f_\mathrm{spa}(h) + f_\mathrm{phy}(h)$.
% 
% * The splitting in $F$ and $G$ is done to be in the setting of primal-dual
% algorithm, see [2]. (The primal problem can be reformulated as primal-dual
% problem which leads to the primal-dual algorithm.)
%
% * Other identifications of the functionals $F$ and $G$ are possible, but in
% public version of this package we only provide this one.
% (This identification is set in <setInvType.html> with parameter
% |seti.invNo = 6|.)
% 
%
% *The primal-dual algorithm*, for convex problems, see Algorithm 1 in [2]:
%
% * Note that $F^\ast$ is the Fenchel conjugate of $F$.
% * The x in the algorithm corresponds to h above.
%
% * Initialization:
%
% Choose
% primal step size $\tau > 0$,
% dual step size $\sigma > 0$,
% initial vectors $(x^0,y^0) \in X \times Y$
% (e.g. $x^0 = 0$, $y^0 = 0$),
% and set $\bar{x}^0 = x^0$.
%
% * Iterations ($n \geq 0$): Update $x^n, y^n, \bar{x}^n$ as follows:
%
% 1. $\quad y^{n+1} = (I + \sigma \partial F^\ast)^{-1}
%  (y^n + \sigma K \bar{x}^n)$,
%
% 2. $\quad x^{n+1} = (I + \tau \partial G)^{-1}
%  (x^n - \tau K^\ast y^{n+1})$,
%
% 3. $\quad \bar{x}^{n+1} = 2x^{n+1} - x^n.$
%
% The motivation of this algorithm is: 
% 3. is an over-relaxation step, and
% 1.-2. are fixed-point iterations (derived from extremality conditions).
%
%
% *Application of primal-dual algorithm in this package*
%
% To apply the primal-dual algorithm in the context of inverse scattering
% we have to remember that the contrast $q$ is complex, but 
% real vector spaces are indeed crucial for the primal-dual algorithm.
%
% Therefore *transformation operators* are introduced in 
% <setIdImagReal.html>, see also Section 4.5 in [1]:
%
% * $\texttt{seti.S}
%  = T_{\bf{C} \to \bf{R}^2} : \bf{C} \to \bf{R} \times \bf{R}$, $\quad$
% $T_{\bf{C} \to \bf{R}^2}(x) = (\mathrm{real}(x),\,  \mathrm{imag}(x))$.
%
% * $\texttt{seti.T} 
% = T_{\bf{R} \to \bf{C}} : \bf{R} \times \bf{R} \to \bf{C}$, $\quad$
% $T_{\bf{R} \to \bf{C}} = 
%  y^{\mathrm{real}} + \mathrm{i}\,y^{\mathrm{imag}} 
%  \mathrm{\ where\ } y = (y^{\mathrm{real}},y^{\mathrm{imag}})$.
%
%
% The full iterative reconstruction scheme 
% (with respect to transformation operators)
% is given in Section 4.4 of [1].
%
% A full derivation is given in Setion 4 of [1].
%
% The defintion of functions is done in <setFuncsPda.html> 
% (i.e. discrepancy, penalty terms, components of K and their adjoints,
% ...).
%
%
% *Grid scaling*
%
% Additionally, in the code a grid scaling is respected to offer to compute
% the reconstruction on a coarser grid. Therefore we use
%
% * seti.GU     :   function to scale up the grid,
% * seti.GD     :   function to scale down the grid.
%
% For details see <setGridScale.html>.
%
% Note that grid scaling is not explained explicitly in [1].
%
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
% * [2] Antonin Chambolle and Thomas Pock. A first-order primal-dual algorithm for convex problems with applications to imaging. _Journal of Mathematical Imaging and Vision_, 40(1):120-145, 2011.
%
%
%% See Also
%
% * <start.html>
% * <minPda.html>
% * <minTolIn.html>
% * <minTolOut.html>
% * <setInvType.html>
% * <setIdImagReal.html>
% * <setFuncsPda.html>
% * <setGridScale.html>
%
%
%% Code: function: pda
%
function [hsolCVU,pdaStopInd,FGval,Fval,Gval,relLinDis,disLinInPda,errInPda,minf] = pda(iOut,qCVU,ThetaiOut,seti,dispDepth)
% pda: primal dual algorithm

if ~isequal(size(qCVU),[seti.nROI^seti.dim,1])
    error('size of qCVU does not fit.')
end
% qROI is a complex vector (and not scaled down)

% qCVU = qROI: complex vector upscaled (full) grid

% q: contrast: fixed in pda
% h: update of contrast (changes in pda) (in minPda: new q = q + h)

% Meaning of characters:
% C, R: stored complex or as real values
% M, V: matrix or vector
% U, D: full grid (upscaled, nROI) or small grid (downscaled, nInv)
% e.g.: CMU... is a complex matrix upscaled

% seti.pdaStepsize = 'adaptive'; % fix or adaptive
% is set in checkConsisRec (or in a file inseti...)

GUCV = seti.GU; % input and output is a complex vector (not stored as matrix)
GDCV = seti.GD; % input and output is a complex vector (not stored as matrix)
% iG and G inside GU and GD does matrix <-> vector

% input and output is the complex vector stored in real (RxR instead of C)
GURV = @(x) seti.S(GUCV(seti.T(x))); % RxR -> C, then GU fits, then C -> RxR
GDRV = @(x) seti.S(GDCV(seti.T(x))); % analog; GDRV is currently unused

% qROI is fix in pda algorithm
% store it in different formats... (to use it fast)
% qCVU = qCVU;
qRVU = seti.S(qCVU); % currently unused
qCVD = GDCV(qCVU);
qRVD = seti.S(qCVD); % currently unused

% qROI not down-scaled in this file
% (and not before; qROI is a fixed vector in this file)

% x is down-scaled
% final h is upscaled

% M1a, M1b only in case of seti.useWavelet == 1

% hs is real x real; h is corresponding complex
% inside pda: hs is used, but output is complex h

% in case of using pda: checkConsistency.m set
%seti.pNorm = 2; % then WHS-Norm is used(!) (important)
%seti.qNorm = 1;
% if you change pNorm or qNorm you have to change proximal mappings etc.

% inside pda h is real: hs = [hr; hi] = [real(h); imag(h)]
% outside: h is complex (so h = hz)

%[m,~,~] = size(A); 1:m/2 is real part and m/2+1:end is imag part

% notation:
% qROI is complex(!)
% qr = real(qROI);
% qi = imag(qROI);
% q = [qr;qi]; % in 2D: in R^{2 nROI^2}; in 3D: in R^{2 nROI^3}
% analog: h

% in this file: q = qROI is qROI and NOT W(qROI)

speed = 0;
stopPda = 0; % 0 or 1: break pda if changes are below tolerance

nPda = seti.pdaN; % number of iterations

reMin = seti.physBounds(1);
reMax = seti.physBounds(2);
imMin = seti.physBounds(3);
imMax = seti.physBounds(4);

% first-order primal-dual algorithm (pda)
% for convex problems (paper: Chambolle, Pock, 2011)

% input:
% K: X -> Y: continuous linear operator
% nonlinear primal problem: min_{x \in X} F(Kx) + G(x)
% FC:= F; GC:= G (to have different names to shrinakge...)

% Symbols:
%
% F^\star: convex conjugate of F
% K^\ast: adjoint of K (symbol: *)

% P = seti.pNorm;
% Q = seti.qNorm;
% setting must be: P = 2; Q = 1; % see setInvType

if dispDepth >= 3
    disp('    - Compute auxiliary matrices JA and JB for Jacobian matrix');
end
tic
[JA,JB] = mimo(seti, qCVU, 'jacobian');
tocJac = toc(tic);
if dispDepth >= 3
    fprintf('      Elapsed time of computation of Jacobian matrices is %05.1f min.\n',tocJac/60);
end

%% Code: function: pda: F and G terms
% tF: term F
% tG: term G
clear tF tG;
% tF and tG are set in setInvType
% look there for details, which terms are F and G (or not used)
tF = seti.tF;
tG = seti.tG;
tFG = [tF,tG]; % term FG contains terms of F and G

%% Code: function: pda: initialization

qCVDsize = size(qCVD); % q complex vector downsize
% Consistency check
if ~isequal(qCVDsize,[seti.nInv^seti.dim,1])
    error('pda: qCVDsize does not fit.')
end

yCMsize = size(seti.FmeasDelta); % y size (complex matrix)
if ~isequal(yCMsize,[seti.measNb,seti.incNb])
    error('pda: ysize does not fit.')
end

xnRVD = [zeros(qCVDsize); zeros(qCVDsize)]; %x0 \in X = R^{2n} % real and down
% consistency check
if ~isequal(size(xnRVD),[2*seti.nInv^seti.dim,1])
    error('pda: xnRVD does not fit.')
end

xnRVDsize = size(xnRVD); % currently unused
xmnRVD = xnRVD;

% -- choose K components and initializing y...
% if fclass(i,j) == 1 then you have to choose the fitting K component
% and define yn
% initialize yd0, ys0, ...
% using ysn in case of fs and fsw (because you will not use fs and fsw)
if ismember('fd',tF) || ismember('fdw1',tF) || ismember('fdw2',tF) || ismember('fdw3',tF) % f_d or f_dw1 or f_dw2 in F
    Kcomp.Kd = seti.Kd;
    ydnRM = [zeros(yCMsize); zeros(yCMsize)]; % y_d_0 \in Y = R^{2m x i}
    % y_d_0 \in Y_1=R^(2 measNb x incNb) (data), short Y=R^{2m x i}
    ydRMsize = size(ydnRM);
end
if ismember('fsw',tF) % f_sw in F
    Kcomp.Ksw = seti.Ksw;
    yswn = xnRVD; % y_s_0 \in Y_2 = R^{2n} (in sparsity term)
end
if ismember('fg',tF) % f_g in F
    Kcomp.Kg = seti.Kg;
    % ygn: size of grad(u) that is stored with components as real values (2*...)
    % in 2D: 2*dim x nInv x nInv
    % in 3D: 2*dim x nInv x nInv x nInv
    
    if seti.dim == 2
        ygnRMD = zeros(2*seti.dim,seti.nInv,seti.nInv); % y_g_0 \in Y_g (in gradient sparsity term)
    elseif seti.dim == 3
        ygnRMD = zeros(2*seti.dim,seti.nInv,seti.nInv,seti.nInv);
    end
    %ygn = zeros(2*seti.dim,qsize); % does not work
end
if ismember('fs',tF)
    error('fs in term F (tF) is not implemented.')
end

% KcomponentsStruct can be outside of pda (is needed one time)
[KcompNorm,KvarNames] = KcomponentsStruct(Kcomp,seti);
if iOut == 1 && dispDepth >= 3
    disp(' ')
    fprintf('    - Using K components: ')
    disp(KvarNames')
end

F = Fsum(tF,seti);
G = Gsum(tG,seti);

%DFFq = @(hz) JA*diag(hz)*JB; % for output: FF'(q)[h] (complex)
if strcmp(seti.model,'helmholtz2D') || strcmp(seti.model,'helmholtz3D')
    DFFq = @(xnRVD) JA*diagsparse(seti.GU(seti.T(xnRVD)))*JB;
elseif strcmp(seti.model,'helmholtzHMode2D')
    DFFq = @(xnRVD) JA(:,:,1)*diagsparse(seti.GU(seti.T(xnRVD)))*JB(:,:,1)...
        + JA(:,:,2)*diagsparse(seti.GU(seti.T(xnRVD)))*JB(:,:,2);
else
    %fprintf(strcat('Error in setFuncs - pda not implemented for model ',  seti.model))
    fprintf('Error in pda - pda not implemented for model %s.\n',  seti.model)
end
% with q = qCVU as input (you can use GU(qCVD) too)
% GUCV(T( . )): RVD -> CVU

if dispDepth >= 3
    disp('    - Compute FF(q)');
end
ticFFq = tic;

FFqMeas = mimo(seti, qCVU, 'simo');

% Slower (factor 2) alternative:
% FFqMeas = mimo(seti, qCVU,seti.T(ydnRM)); % expects complex qROI, Fmeas is 0; FFq is complex
% This is alternative is slower, because adjOfDer is computed but not
% needed.

% Notation in general:
% FFqmF  := FFqMeas - FmeasDelta
% ADFFq := FF'(q)^*[FF(q) - FmeasDelta] (Adjoint of derivative)
tocFFq = toc(ticFFq);
if dispDepth >= 3
    fprintf('      Elapsed time of FF(q) computation is %05.1f min.\n',tocFFq/60);
end

%% Code: function: pda: Algorithm 1

% ynp = y_{n+1}
% xnp = x_{n+1}
% xmnp = \={x}_{n+1} (x with makron)

% yn = y_n
% xn = x_n
% xmn = \={x}_n (x with Makron)

% x     : primal variable
% tau   : primal stepsize (tau > 0)
% y     : dual variable
% sigma : dual stepsize (sigma > 0)

% exact: ||K|| = max_{||h||_1 \neq 0} ||Kh||_HS / ||h||_1
% implementation with random Vektors h...

if dispDepth >= 3
    disp('    - L operator norm: numerical approximation')
end
ticOpNorm = tic;
L = opNormNum(xnRVD,JA,JB,KcompNorm,seti,dispDepth);
tocOpNorm = toc(ticOpNorm);
if dispDepth >= 3
    fprintf('      Elapsed time of operator norm L is %05.1f min.\n',tocOpNorm/60);
end
clear ticOpNorm tocOpNorm;
if dispDepth >= 3
    fprintf('      L = ||K|| = %g\n',L);
end

theta = 1; % pda algorithm theta = 1 such that the proof works (Th. 1 in [2].)
vartheta = seti.vartheta; % vartheta \in (0,1) in choosingStepsizes!!! (not theta from PDA...)

% pda algorithm

%% Code: function: pda: control results of F(K(h)) and G(h) (no influence on calculation)
Fval = zeros(nPda,1); % store values of FC(Kh)
Gval = zeros(nPda,1); % store values of GC(h)
FGval = zeros(nPda,1);
disLinInPda = ones(nPda,1); % linearized relative discrepancy
% disLinInPda = \|FF'(q)[h]+FF(q)-F_\meas^\delta\|_{WS,pNorm} / \|y^\delta\|_{WS,pNorm}
relLinDis = zeros(nPda,1); % relLinDis = disLin / dis;
errInPda = zeros(nPda,1); % relative error in pda

fd = zeros(nPda,1);
fs = zeros(nPda,1);
fg = zeros(nPda,1);
fp = zeros(nPda,1);

% Define variables... will be set in pdaChoosingStepsizes
xnUp = GURV(xnRVD);
xnpRVDstore = xnUp; % choosing needs upscaled
xnRVDstore = xnUp;
tau = 0;
sigma = 0;
tauVal = zeros(nPda,1); % vector to store tau values
sigmaVal = zeros(nPda,1); % vector to store sigma values

%% Code: function: pda: pda
for iPda = 1:nPda
    
    [tau,sigma] = pdaChoosingStepsizes(seti.pdaStepsize,L,iPda,xnpRVDstore,xnRVDstore,KcompNorm,tau,sigma,vartheta,seti,JA,JB);

    tauVal(iPda) = tau; % just to store...
    sigmaVal(iPda) = sigma; % just to store
    
    %----- pda algorithm steps -----

    %-- definition
    %wyi = @(yin,Ki) yin + sigma*GDRV(Ki(GURV(xmnRVD))); % down scaled (not for wyd)
    wyi = @(yin,Ki) yin + sigma*Ki(xmnRVD);
    % i is replaced by d, s or g (d and g extra cases...)
    wyd = @(ydn,Kd) ydn + sigma*Kd(xmnRVD,JA,JB);
    wyg = @(ygn,Kg) ygn + sigma*Kg(xmnRVD);
    wxF = @(KAdjy) xnRVD-tau*KAdjy; % w_x function handle
    %wx = wxF(KAdjy); % how to call it
    
    %-- compute KAdjySum (prepare pda step 2)
    KAdjySum = 0;

    %-- pda step 1: prox_F^+ (resolvent to F^+)
    % and prepare pda step 2 (compute KAdjySum)
    % 1) fd etc. in F (a discrepancy term should be in F)
    
    wydVal = wyd(ydnRM,seti.Kd); % w_{y_d}
    if ismember('fd',tF) % using Kd in cases of fd
        ydnp = proxFdPlus(sigma,seti.vd(FFqMeas),wydVal); % prox_{sigma F_d^+}(w_y_d)
        %fprintf('proxFdPlus: sigma = %g | max(vd) = %g | max(wyd) = %g\n',sigma,max(max(abs(seti.vd(FFq)))),max(max(abs(wydVal))));
    elseif ismember('fdw1',tF) % fdw1 in F (case seti.p = 1)
        ydnp = proxFdw1Plus(sigma,seti.vd(FFqMeas),wydVal,seti); % prox_{sigma F_d^+}(w_y_d)
    elseif ismember('fdw2',tF) % fdw2 in F (case seti.p == 2)
        ydnp = proxFdPlus(sigma,seti.vd(FFqMeas),wydVal); % prox_{sigma F_d^+}(w_y_d)
        % proxFdw2Plus = proxFdPlus because fdw2 = fd
    elseif ismember('fdw3',tF) % fdw3 in F (case seti.pNorm ~= 2)
        ydnp = proxFdw3Plus(sigma,seti.vd(FFqMeas),wydVal,seti); 
        % proxFdw2Plus = proxFdPlus because fdw2 = fd
    else
        error('No discrepancy term choosen... Does not make sense.')
    end
    KdAdjyd = seti.KdAdj(ydnp,JA,JB); % preparation of pda step 2

    if dispDepth >= 5 % additional information
        fprintf('   norm(KdAdjyd) = %g\n',norm(KdAdjyd));
        fprintf('   max(KdAdjyd) = %g\n',max(KdAdjyd));
        fprintf('   min(KdAdjyd) = %g\n',min(KdAdjyd));
    end
    
    % -- test start
    %seti.KdAdj = @(yd,JA,JB) seti.S(seti.GD(ADFFqFast(seti.T(yd),JA,JB,seti)));
    
    if 0
    if iPda == 10
        disp('iPda = 10');
        
        figure(101);
        imagesc(real(seti.T(ydnp))); colorbar; axis xy;
        figure(102);
        imagesc(imag(seti.T(ydnp))); colorbar; axis xy;
        
        figure(103);
        imagesc(real(seti.T(ydnp))+imag(seti.T(ydnp))); colorbar; axis xy;
        
        A = ADFFqFast(seti.T(ydnp),JA,JB,seti);
        Amat = seti.G(A);
        figure(104);
        imagesc(real(Amat)); colorbar; axis xy;
        figure(105);
        imagesc(imag(Amat)); colorbar; axis xy;
        
        B = seti.GD(A);

        C = seti.S(B);

        figure(106);
        imagesc(seti.G(seti.R(xnRVD))); colorbar; axis xy;
        figure(107);
        imagesc(seti.G(seti.I(xnRVD))); colorbar; axis xy;

        error('stop');
    end
    end
    
    % -- test end
    
    KAdjySum = KAdjySum + KdAdjyd; % preparation of pda step 2

    % 2) fs etc. in F
    if ismember('fs',tF) % fs in F
        if dispDepth >= 4
            disp('case fs in F not yet.')
        end
    elseif ismember('fsw',tF) % fsw in F
        wysw = wyi(yswn,seti.Ksw); % w_{y_sw}
        yswnp = proxFswPlus(sigma,seti.vsw(qCVU),wysw,seti); % prox_{sigma F_s^+}(w_y_s)
        % ysnp is down scaled
        %KswAdjys = seti.KswAdj(GURV(yswnp)); % preparation of pda step 2
        KswAdjys = seti.KswAdj(yswnp); % preparation of pda step 2
        KAdjySum = KAdjySum + KswAdjys; % preparation of pda step 2
        clear wysw KswAdjys;
    else
        if iOut == 1 && iPda == 1
            if dispDepth >= 4
                disp('    - Note: No term like fs (sparsity) is used in F.')
            end
        end
    end


    % 3) fg etc. in F
    if ismember('fg',tF) % fg in F
        wygVal = wyg(ygnRMD,seti.Kg); % w_{y_g}
        ygnp = proxFgPlus(sigma,seti.vg(qCVU),wygVal,seti); % prox_{sigma F_g^+}(w_y_g)
        %fprintf('proxFgPlus: sigma = %g | max(vg) = %g | max(wyg) = %g\n',sigma,max(max(max(abs(seti.vg(qCVU))))),max(max(max(abs(wygVal)))));
        KgAdjyg = seti.KgAdj(ygnp); % preparation of pda step 2
        
        if dispDepth >= 5 % additional information
            fprintf('   norm(KgAdjyg) = %g\n',norm(KgAdjyg));
            fprintf('   max(KgAdjyg) = %g\n',max(KgAdjyg));
            fprintf('   min(KgAdjyg) = %g\n',min(KgAdjyg));
        end

        KAdjySum = KAdjySum + KgAdjyg; % preparation of pda step 2
        %clear wygVal KgAdjyg;
    end
    
    % 4) fp etc. in F
    if ismember('fp',tF)
        if dispDepth >= 4
            disp('fp in F not implemented. Maybe not useful.')
        end
    end

    %-- pda step 2: prox_G(w_x) (Resolvente zu G)
    wx = wxF(KAdjySum); % grid down and real

    % order: dsgp, e.g. proxGsp, proxGp
    % fitting to fclass = [<d> <s> <g> <p>]
    % f_d in G does not make sense

    %fG = (fclass == 2);
    if ismember('fd',tG) || ismember('fdw1',tG) || ismember('fdw2',tG)
        error('f_d in Term G does not make sense.')
    elseif isequal({'fs'},tG)
        xnpRVD = proxGs(tau,seti.alpha,seti.vs(qCVU),wx); % prox_{\tau G_s}(w_x)
        % to do: implement proxGs
    elseif isequal(sort({'fs','fp'}),sort(tG))
        vsp = seti.vs(qCVU); % seti.vs = seti.vp in setiFuncs
        % qCVU does not change...
        % you can compute it at the beginning of pda
        % the same with vp, vsp...
        xnpRVD = proxGsp(tau*seti.alpha*seti.dVinv,vsp,wx,seti);

        if dispDepth >= 5 % additional information
            fprintf('   max(wx) = %g\n',max(wx));
            fprintf('   max(xnpRVD) = %g\n',max(xnpRVD));
            fprintf('   kappa = %g\n',tau*seti.alpha*seti.dVinv);
        end
        %fprintf('proxGsp: tau = %g | alpha*dV = %g | max(wx) = %g | max(vsp) = %g\n',tau,seti.alpha*seti.dVinv,max(max(abs(wx))), max(abs(vsp)));
        %xnp = proxGsp(tau,alpha,vsp,wx); % prox_{\tau G_sp}(w_x)
    elseif isequal(sort({'fs','fp2'}),sort(tG))
        vsp = seti.vs(qCVU);
        xnpRVD = proxGsp2(tau*seti.alpha*seti.dVinv,vsp,wx,seti);

        if dispDepth >= 5 % additional information
            fprintf('   max(wx) = %g\n',max(wx));
            fprintf('   max(xnpRVD) = %g\n',max(xnpRVD));
            fprintf('   kappa = %g\n',tau*seti.alpha*seti.dVinv);
        end
    elseif isequal({'fp'},tG)
        % to do: rewrite seti.vsw, so that the arg is grid down
        vp = seti.vs(qCVU);
        xnpRVD = proxGp(reMin,reMax,imMin,imMax,vp,wx,seti); % prox_{\tau G_p}(w_x)
    else
        disp('The wanted prox_G is not implemented yet.')
        error('stop');
    end
    % xnp: grid down and real

    %-- pda step 3:
    xmnpRVD = xnpRVD+theta.*(xnpRVD-xnRVD); % xmnp: grid down and real

    %----- end of: pda algorithm steps -----

    
    % store to use in adaptive stepsize (is needed upscaled)
    xnRVDstore = xnRVD;
    xnpRVDstore = xnpRVD;

    %-- prepare new iteration
    
    if exist('ydnp','var') && isnumeric(ydnp)
        ydnRM = ydnp;
    end
    if exist('yswnp','var') && isnumeric(yswnp)
        yswn = yswnp;
    end
    if exist('ygnp','var') && isnumeric(ygnp)
        ygnRMD = ygnp;
    end
    if exist('ypnp','var') && isnumeric(ypnp)
        ypn = ypnp;
    end

    xnRVD = xnpRVD;
    xmnRVD = xmnpRVD;
    
    if speed == 0

        % xnRVD: solution (real vector grid down)
        
        %-- store values of functions F, G

        DFFqh = DFFq(xnRVD);
        Fval(iPda) = F(DFFqh,FFqMeas,qCVU,xnRVD);
        Gval(iPda) = G(qCVU,xnRVD);
        FGval(iPda) = Fval(iPda) + Gval(iPda);

        %relLinDisc = linDiscr(DFFqh,FFq)/linDiscr(0*DFFqh,FFq); % altern. for abs lin. dis.
        %relLinDisc = linDiscr(DFFqh,FFq)/normws(FFq-seti.FmeasDelta,seti);
        
        disLinInPda(iPda) = normws(seti.S(DFFqh)+seti.S(FFqMeas-seti.FmeasDelta),seti)/normws(seti.FmeasDelta,seti);
       
        % dis of last iOut
        if iOut == seti.iOutIni+1
            dis = seti.disIni;
        else
            dis = seti.dis(iOut-1);
        end
        relLinDis(iPda) = disLinInPda(iPda)/dis;
        clear dis;
        hsolCVUcurrent = seti.GU(seti.T(xnRVD));
        errInPda(iPda) = norm(seti.qROIexact-(qCVU+hsolCVUcurrent),2)/norm(seti.qROIexact,2); % relative error in 2-Norm
       
        %-- compute discrepancy and penalty terms

        if dispDepth >= 4
            fprintf('    - i = %03d | stepsizes: tau = %1.2g, sigma = %1.2g\n',iPda,tau,sigma)
            fprintf('               F = %1.2g | G = %1.2g | F + G = %1.2g \n',Fval(iPda),Gval(iPda),FGval(iPda))
        end
        if ismember('fd',tFG)
            fd(iPda) = seti.fd(DFFqh,FFqMeas);
            if dispDepth >= 4
                fprintf('               f_d = %1.2g \n',fd(iPda))
            end
        end
        if ismember('fs',tFG)
            fs(iPda) = seti.fs(qCVU,xnRVD);
            if dispDepth >= 4
                fprintf('               f_s = %1.2g\n',fs(iPda))
            end
        end
        if ismember('fg',tFG)
            fg(iPda) = seti.fg(qCVU,xnRVD);
            if dispDepth >= 4
                fprintf('               f_g = %1.2g\n',fg(iPda))
            end
        end
        if ismember('fp',tFG)
            fp(iPda) = seti.fp(qCVU,xnRVD);
            if dispDepth >= 4
                fprintf('               f_p = %1.2g\n',fp(iPda))
            end
        end
        if dispDepth >= 4
            fprintf('               disLinInPda = %1.2g | RelLinDis = %1.2g \n',disLinInPda(iPda),relLinDis(iPda))
            fprintf('               errInPda = %1.2g \n',errInPda(iPda))
        end
       
        if iPda > 0 && seti.useTolIn == 1 && relLinDis(iPda) < ThetaiOut
            if dispDepth >= 3
                disp('    - Break pda because inner tolerance principle.')
            end
            break;
        end
        
        if stopPda == 1 && iPda >= 3 && norm(FGval(iPda)-FGval(iPda-2)) + norm(FGval(iPda)-FGval(iPda-1)) < 1E-6*norm(FGval(iPda))
            if dispDepth >= 3
                disp('    - Break pda because changes are below tolerance.')
            end
            break;
        end

    end
    
    % plots inside pda not available in public version
    if seti.plotFreqiPda ~= 0 && (seti.plotFreqiPda == 1 || floor(iPda/seti.plotFreqiPda) == iPda/seti.plotFreqiPda)
        plotInsidePda;
    end

    % --

    % final:
    % x = xnp: min_x F(Kx) + G(x); here: x = hs;
    
end

pdaStopInd = iPda;

if speed == 1
    % xnRVD is the solution h (but grid down and real)
    xnCVD = seti.T(xnRVD);
    DFFqh = DFFq(xnCVD);
    Fval(pdaStopInd-1) = F(DFFqh,FFqMeas,qCVU,xnCVD);
    Gval(pdaStopInd-1) = G(qCVU,xnCVD);
    FGval(pdaStopInd-1) = Fval(iPda) + Gval(iPda);
end

if speed == 0
    % parts of functional f to minimize
    minf.fd = fd;
    minf.fs = fs;
    minf.fg = fg;
    minf.fp = fp;
end

% output solution: xnRVD -> xnCVU
hsolCVU = seti.GU(seti.T(xnRVD)); % Update (q := q+h) (q,h complex)

end

%% Code: subfunction to K components norm
%
% *KcomponentsStruct*
%
function [KcompNorm,KvarNames] = KcomponentsStruct(Kcomp,seti)
% like Kcomponents but with input Kcomp as struct
%n = length(fieldnames(Kcomp)); % number of K_1, K_2, ..., K_n
fields = fieldnames(Kcomp); % e.g. 'Kd' 'Kg' 'Kp'; varNames(1) = 'Kd'
% if Kd... then normws2, otherwise normroi2
KcompNorm = @(xnRVD,JA,JB) 0;
% i = 1: must be Kd, so with normws2 norm
i = 1;
KcompNorm = @(xnRVD,JA,JB) KcompNorm(xnRVD,JA,JB) + normws2(Kcomp.(fields{i})(xnRVD,JA,JB),seti)^2; 
for i = 2:numel(fields)
    %Kcomp.(fields{i})
    if strcmp('Kg',fields(i))
        KcompNorm = @(xnRVD,JA,JB) KcompNorm(xnRVD,JA,JB) + normTVinv2(Kcomp.(fields{i})(xnRVD),seti)^2;
    else
        KcompNorm = @(xnRVD,JA,JB) KcompNorm(xnRVD,JA,JB) + norminv2(Kcomp.(fields{i})(xnRVD),seti)^2;
    end
    %fprintf('Kcomp inputname: %s\n',inputname(i))
    %KvarNames{i} = fields{i};
end
KcompNorm = @(xnRVD,JA,JB) sqrt(KcompNorm(xnRVD,JA,JB));
% KcompNorm = \sqrt{ \|\cdot\|_WHS^2 + \|\cdot\|_\ROI^2 + ... }
KvarNames = fields;
end

%% Code: subfunctions to F and G (sum)
%
% *Fsum*
%
function F = Fsum(tF,seti)
% e.g. F = f_d + f_s, when tF stores fd and fs
% call F = Fsum(tF,seti) and later F(DFFqh,FFq,qCVU,xnRVD)
%old q = qCVU, old hz = xnRVD
Fsum = @(DFFqh,FFq,qCVU,xnRVD) 0;
for i = 1:length(tF)
    
    fj = seti.(tF{i});
    if ismember(tF(i),{'fd','fdw1','fdw2','fdw3'}) % f_d, f_dw etc. does need input arguments DFFqh and FFq
        Fsum = @(DFFqh,FFq,qCVU,xnRVD) Fsum(DFFqh,FFq,qCVU,xnRVD) + fj(DFFqh,FFq);
    else % other like f_s, f_p does need input arguments q and hz
        Fsum = @(DFFqh,FFq,qCVU,xnRVD) Fsum(DFFqh,FFq,qCVU,xnRVD) + fj(qCVU,xnRVD);
    end
    clear fj;
end
%if strcmp('tF',inputname(1)) % input tF, so you want to get F
F = @(DFFqh,FFq,qCVU,xnRVD) Fsum(DFFqh,FFq,qCVU,xnRVD);
clear i Fsum;
end

%%
% *Gsum*
%
function G = Gsum(tG,seti)
% similar to Fsum, but only 2 arguments
% call G = Gsum(tG,seti)
Gsum = @(q,hz) 0;
for i = 1:length(tG)
    fj = seti.(tG{i});
    Gsum = @(q,hz) Gsum(q,hz) + fj(q,hz);
    clear fj;
end
G = @(q,hz) Gsum(q,hz);
clear i Gsum;
end

%% Code subfunctions to proximal mappings
%
%%
% *proxFdPlus*
%
function res = proxFdPlus(sigma,vd,wyd) % d: discrepancy
% f_d(x) = 1/2 ||S(FF'(q)[x]) + vd||_WS,2^2, here with vd = S(FFq-y^\delta)
% res = y_d^(n+1) = \prox_{\sigma F_d^+} (w_y_d)
res = (wyd+sigma*vd)/(1+sigma);
end

%%
% *proxFdw1Plus*
%
function res = proxFdw1Plus(sigma,vd,wyd,seti)
%f_dw(x) = 1*||S(FF'(q)[x]) + vd||WS,2^1, here with vd = S(FF(q)-y^\delta)
% res = y_dw1^(n+1) = \prox_{\sigma F_dw1^+} (w_y_dw1)
% seti in normws is seti.pNorm, here seti.pNorm must be 2
% (in pda with invOps=3,4,5,6 always seti.pNorm = 2 and seti.qNorm = 1)
res = (wyd+sigma*vd)/max(1,normws(wyd+sigma*vd,seti));
end

%%
% *proxFdw3Plus*
%
function res = proxFdw3Plus(sigma,vd,wyd,seti)
%f_dw3(x) = (1/p)*||S(FF'(q)[x]) + vd||_p^p, here with vd = S(FF(q)-y^\delta)
% res = y_dw1^(n+1) = \prox_{\sigma F_dw3^+} (w_y_dw1)

if seti.pNorm > 1
    % requires seti.qNorm>1
    pPrime = seti.pNorm/(seti.pNorm-1);
    res = shrinkFuncComp(real(wyd+sigma*vd),sigma,pPrime,'newton')...
        +1i*shrinkFuncComp(imag(wyd+sigma*vd),sigma,pPrime,'newton'); 
        
elseif seti.pNorm==1
    res = max((abs(real(wyd+sigma*vd))-sigma),0).*sign(real(wyd+sigma*vd)) ...
        +1i*max((abs(imag(wyd+sigma*vd))-sigma),0).*sign(imag(wyd+sigma*vd)); 
else
    disp('seti.qNorm is less than one in pda.m (?!)');
end

end

%%
% *proxFswPlus*
%
function res = proxFswPlus(sigma,vs,wys,seti) % s: sparse
% f_s(x) = ||alpha*x + v_s||_ROI,1
% factor alpha must be inside v_s
% res = y_s^(n+1) = \prox{\sigma F_s^+} (w_y_s)
rho = seti.dVinv; % This factor was corrected from 1 to seti.dVinv on 20181213.
res = intProj(-rho,rho,wys+sigma*vs);
end

%%
% *proxFgplus*
%
function res = proxFgPlus(sigma,vg,wyg,seti) % g gradient
% f_g(x) = ||beta*\grad(x) + v_g||_ROI,1
% factor beta must be inside v_g
% res = y_g^(n+1) = \prox{\sigma F_g^+} (w_y_g)
% --
%res = intProj(-1,1,wyg+sigma*vg); % used until 20160406, and [20160412,20160414] % wrong...
%new: 20160406--20160412: was wrong...:
%t = wyg+sigma*vg;
%res = t./max(1,abs(t)/seti.dV);
% --
% new...

z = wyg + sigma*vg;
%prepare matrix of norm of tuples...
d = size(z,1)/2;
nInv = size(z,2);
sizez = size(z);
zreshape = reshape(z,[2*d nInv^d]); % for 2D and 3D
% compute norms
tupelnorms = transpose(sqrt(sum(abs(zreshape).^2,1)));
%tupelnormsReshape = reshape(tupelnorms,sizey(2:end)); % for 2D and 3D
tupelnorms = reshape(tupelnorms,[1 sizez(2:end)]); % second reshape to add one dimension
tupelnorms = repmat(tupelnorms,[2*d 1 1]); % repmat for element-wise use...

% tupelnormsMatElWise has the the same dimension as z and can used component-wise
z = wyg+sigma*vg;
res = z./max(1,tupelnorms);
end

%%
% *proxGsp*
%
function res = proxGsp(kappa,vsp,wx,seti) % s: sparse AND p: physical
% G(x) = f_s(x) + f_p(x)
%      = \alpha ||x + v_s||_1 + \delta_[a,b](Re(x + v_p)) + \delta_[c,d](Im(x + v_p))
% important: v_s == v_p
% res = x^(n+1) = \prox_{\tau G_sp} (w_x)
% (shrinkage with bounds a and b)
%res = [shrinkage(wx+vsp,tau*alpha,a,b); shrinkage(wx+vsp,tau*alpha,c,d)]-1/seti.alpha*seti.S(vsp);

if isfield(seti,'omegaW1p')
    if strcmp(seti.wavIsom,'W1p') 
        res = [seti.shrkRe(seti.R(wx+vsp),kappa*seti.omegaW1p); seti.shrkIm(seti.I(wx+vsp),kappa*seti.omegaW1p)]-vsp;
    else
        error('Error in pda.m: seti.wavIsom not implemented.') 
    end
else
    res = -vsp+[seti.shrkRe(seti.R(wx+vsp),kappa); seti.shrkIm(seti.I(wx+vsp),kappa)];
end

end

%%
% --------------------------------------------------------------------
% *proxGsp2* (fs and fp2 in G)
%
function res = proxGsp2(kappa,vsp,wx,seti) % s: sparse AND p: physical (p2 distinguishes obstacle and background)
% This function is similar to proxGsp. The difference is the splitting of
% the contrast in the obstacle and the background.
%
% G(x) = f_s(x) + f_p2(x)
%      = \alpha ||x + v_s||_1 + \delta_[a,b](Re(obs(x + v_p))) + \delta_[c,d](Im(obs(x + v_p)))
%                             + \delta_[0,0](Re(back(x + v_p))) + \delta_[0,0](Im(back(x + v_p)))
% Notation: obs does mean seti.qObs, and back does mean seti.qBack.
% In mathematical notation we chose: obs(x) = \check{x} and back(x) = \mathring{x}.
% important: v_s == v_p
% res = x^(n+1) = \prox_{\tau G_sp2} (w_x)

% long:
% wvBack = seti.qBack(seti.R(wx+vsp));
% wvBack = zeros(size(wvBack)); % physical bounds of background are zero, i.e. wvBack is zero
% short:
wvBack = zeros(size(seti.nROI^seti.dim,1)); % physical bounds of background are zero, i.e. wvBack is zero

res = -vsp+[seti.qMerge(seti.shrkRe(seti.qObs(seti.R(wx+vsp)),kappa),wvBack); seti.qMerge(seti.shrkIm(seti.qObs(seti.I(wx+vsp)),kappa),wvBack)];

end
% --------------------------------------------------------------------
%%
% *proxGp*
%
function res = proxGp(a,b,c,d,vp,wx,seti) % p: physical... Potentialtopf-Funktion
% f_p(x) = \delta_[a,b](R(wx + v_p)) + \delta_[0,c](I(wx+v_p))
% res = x^(n+1) = \prox_{\tau G_p} (w_x)
res = [intProj(a,b,seti.R(wx+vp)); intProj(c,d,seti.I(wx+vp))]-vp;
end

