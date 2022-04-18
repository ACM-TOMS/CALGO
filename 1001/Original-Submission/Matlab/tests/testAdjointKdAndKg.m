%% Test adjoints of Kd and Kg
% Test of the adjoints of Kd and Kg.
%
%% Syntax
%
%  testAdjointKdAndKg
%
%% Description
%
% |testAdjointKdAndKg| tests numerical whether 
% seti.KdAdj is the adjoint of seti.Kd
% and seti.KgAdj is the adjoint if seti.Kg.
%
% Further we test, whether the norm |normTVinv2| is induced by |innergrad|.
%
% This test should run in environment runtests, see <runtests.html> via
%
%   runtests('a')
%
%% More About
%
% Kd, Kg and their adjoints are explained in <setFuncsPda.html> 
% (seti.Kd, seti.KdAdj, seti.Kg, and seti.KgAdj).
%
% For the inner products
% $\langle \cdot, \cdot \rangle_{\mathrm{dis},\bf{R}}$, 
% $\langle \cdot, \cdot \rangle_{\mathrm{roi},\bf{R}}$, and
% $\langle \cdot, \cdot \rangle_{\mathrm{tv},\bf{R}}$
% see [1, Sec. 4.5].
%
% To test the adjoints we compare:
%
% * $\langle \texttt{Kd}(\texttt{xnRVD}), \texttt{yd}\rangle_{\mathrm{dis},\bf{R}}$
% with $\langle\texttt{xnRVD}, \texttt{Kd}^\ast(\texttt{yd})\rangle_{\mathrm{roi},\bf{R}}$
% * $\langle\texttt{Kg}(\texttt{xnRVD}), \texttt{yg}\rangle_{\mathrm{tv},\bf{R}}$
% with $\langle\texttt{xnRVD}, \texttt{Kg}^\ast(\texttt{yg})\rangle_{\mathrm{roi},\bf{R}}$.
%
% Of course they should be numerical identical.
%
% * This function could be shorter and faster:
%   the functions setData and setRecon does set more fields 
%   than necessary for the test.
% * This function stores the predefined contrast and measurement set-up 
%   in the folder "output". 
%   A flag for <setInput.html> and <setGeomSim.html> to prohibit plotting
%   would be useful.
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <runtests.html>
% * <setFuncsPda.html>
%
%% Code

%%
% *Set it only in case of standalone running this file*

%  clear variables;
%  init;

%%

% Set so that it runs...
seti.invNo = 6;

% set and change to test adjoint...
seti.useWavelet = 0; % or 0... (test both cases)
seti.p = 2;

% set to change...
seti.nCD = 256;

% Load data structure and functions (suppress output by evalc)
disp(' - load preliminary settings with suppressed output')
seti = setGeomSim(seti);
seti = setRecon(seti);
disp(' ')

% Some dimensions for test vectors
n = seti.nROI;
d = seti.dim;
m = seti.measNb;
i = seti.incNb;

% Define test vectors (use dimensions of ROI as vector):
xnRVD = rand(2*n^d,1); % is real (other name: h = xnRVD)
qCVU = rand(n^d,1) + 1i*rand(n^d,1); % is complex (this is qROI = qCVU)

% From pda:
disp('- Jacobian');
tic
[JA,JB] = mimo(seti, qCVU, 'jacobian');
toc

disp('- FF(q)');
tic
FFqMeas = mimo(seti, qCVU, 'simo'); % expects complex qROI
toc

% WHS (i.e. dis,R) and tv,R
yd = rand(size(seti.Kd(xnRVD,JA,JB)));
yg = rand(size(seti.Kg(xnRVD)));

%%
% *Compute scalar products to test adjoints*

% Kd and Kd^\ast
s1 = innerhs(seti.Kd(xnRVD,JA,JB),yd,seti);
s2 = innerroi(xnRVD,seti.KdAdj(yd,JA,JB),seti);
diff = abs(s1-s2)/abs(s1);
fprintf('    diff (rel.) = %g \n',diff);
fprintf('    <Kd(xnRVD), yd>_dis,R = %g+%gi | <xnRVD, Kd^*(yd)>_roi,R = %g+%gi \n',real(s1),imag(s1),real(s2),imag(s2));

% Kg and Kg^\ast
s1 = innergrad(seti.Kg(xnRVD),yg,seti);
s2 = innerroi(xnRVD,seti.KgAdj(yg),seti);
diff = abs(s1-s2)/abs(s1);
fprintf('    diff (rel.) = %g \n',diff);
fprintf('    <Kg(xnRVD), yg>_tv,R = %g+%gi | <xnRVD, Kg^*(yg)>_roi,R = %g+%gi \n',real(s1),imag(s1),real(s2),imag(s2));

%%
% *Test, whether the norm |normTVinv2| is induced by |innergrad|*
%
% Test: Is |normTVinv2(x,seti) = sqrt(abs(innergrad(x,x,seti)))|?

disp(' ')
disp('-- Test: is normTVinv2 induced by innergrad? Yes: ')

s1 = sqrt(abs(innergrad(yg,yg,seti)));
s2 = normTVinv2(yg,seti);
diff = abs(s1-s2)/abs(s1);
fprintf('    diff (rel.) = %g \n',diff);
fprintf('    ||yg||_tv,R,2 = %g+%gi | sqrt(abs(<yg,yg>_tv,R)) = %g+%gi \n',real(s1),imag(s1),real(s2),imag(s2));
