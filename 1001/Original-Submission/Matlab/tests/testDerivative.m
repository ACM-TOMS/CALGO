%% testDerivative
% Test the Frechet derivative and its adjoint of the 
% contrast-to-measurement operator F.
%
%% Syntax
%
%   testDerivative
%
%% Description
% 
% |testDerivative| tests the Frechet derivative and its adjoint of the 
% contrast-to-measurement operator $\mathcal{F}$ by plotting results of 
% the Frechet condition. Further the adjoint of the Frechet derivative is
% tested.
%
% This test should run in environment runtests, see <runtests.html> via
%
%   runtests('d')
%
%% Input Arguments
%
% * seti.dim    :   dimension of the problem (2 or 3)
%
%% Output Arguments
%
% *Figures*, see <runtests.html>.
%
%% More About
%
% * Note the divergent results for single and double precision computation.
%  (Single precision is not computed in case of 3D to save time.)
%
% *Frechet derivative*
%
% The forward operator $\mathcal{F}$ is Frechet differentiable at $q$ if
%
% $\lim_{t \rightarrow 0} 
% \frac{\|\mathcal{F}(q + th) - \mathcal{F}(q) - \mathcal{F}'(q)[th]\|_W}
% {\|th\|_V} = 0$, where
%
% * W = |normws2|, i.e. the norm which is induced by |innerhs|, 
%   i.e. $\langle A,B \rangle_\mathrm{dis}$,
% * V = |normroi2|, i.e. the norm which is induced by |innerroi|, 
%   i.e. $\langle x, y \rangle_\mathrm{roi}$.
%
% We plot the results for several small values of $t$.
%
% Note that $\mathcal{F}'(q)[\cdot]$ is linear by definition, such that
% $\mathcal{F}'(q)[th] = t \mathcal{F}'(q)[h]$, which is implemented.
%
% *Test adjoint of Frechet derivative by inner products*
%
% Compare
% $\langle \mathcal{F}'(q)^ast[\mathcal{F}(q)-y], g \rangle_\mathrm{roi}$
% and 
% $\langle \mathcal{F}(q)-y, \mathcal{F}'(q)[g] \rangle_\mathrm{dis}$ 
% for a complex test vector $g$.
%
% *Test adjoint of the Frechet derivative by direct computation*
%
% Compare the adjoint of the Frechet derivative |ADFFq| computed by <mimo.html> 
% with the result, directly computed by A,B decomposition. 
% See <ADFFqFast.html> for the method.
%
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <runtests.html>
% * <ADFFqFast.html>
%
%% Code

%%
% Only in case of standalone:

%  init;

%%
% *Set parameters --------------------------------------------------------*

% Set dimension of the problem if not already set
if ~isfield(seti,'dim')
    seti.dim = 2;
end

useSmoothTest = 1; % use a smooth test vector instead of rand... (better convergence to 0)

seti.fileSuffix = sprintf('_testDerivative%dD',seti.dim);

seti.rCD = 0.2; % i.e. in notation of [1] R = rCD/2 = 0.1
if seti.dim == 2
    seti.nCD = 256; 
else
    seti.nCD = 64;  % mimo with adjOfDer takes 7 mins for 3D in case of nCD = 256.
                    % Therefore use nCD = 64; later it will be a little bit higher.
end

% -- setGridScale
seti.gscale = 0;
seti.nCDinv = 40;

% -- setKernel
seti.k = 250;
seti.model = 'helmholtz';

% -- setContrast
if seti.dim == 2
    seti.contrast = 'twoCornersOneBall2D'; % this is an interesting case because single and double precision differs
else
    seti.contrast = 'referenceBall'; %auto 2D or 3D
end

% -- setMeas (measurements)
seti.incType = 'pointSource';
seti.measType = 'nearField';
seti.incNb = 35;
seti.measNb = 35;
seti.radSrc = 5;
seti.radMeas = 5;

%%
% *setData ---------------------------------------------------------------*

disp('# setGeomSim')
seti = setGeomSim(seti); % set geometry and simulation

%%
% *make F(q) -------------------------------------------------------------*
y = zeros(seti.measNb, seti.incNb);
disp('Compute FF(q) and ADFFq')
tic
[FFqMeas,ADFFq] = mimo(seti,seti.qROIexact,y,'adjOfDer'); % FFqMeas = FF(q) and ADFFq = DG = FF'(q)*[FF(q)-y]
toc

%%
% *Make Jacobi Matrices s.t. FF'(q)h = JA*diag(h)*JB ---------------------*
disp('make Jacobi Matrices')
tic
[JA,JB] = mimo(seti,seti.qROIexact,'jacobian');
toc

%%
% *Check plots for bugfixing (currently unused)*

% To plot figures on CD...:
E = @(x) extendROItoCD(x,seti.ROImask);

if 0
    figure(101); % OK
    imagesc(real(E(seti.qROIexact))); colorbar; grid on;
    axis xy;
    title('changestructure')
end

if 0
    figure(102); % OK
    gridax = seti.grid(2,1:seti.nCD);
    imagesc(gridax,gridax,real(E(seti.incField(:,1)))); colorbar;
    axis xy;
    title('changestructure')
end

if 0
    figure(103); % OK
    gridax = seti.grid(2,1:seti.nCD);
    
    xplot = seti.k^2.*seti.measKer(4,:);
    
    imagesc(gridax,gridax,real(E(xplot))); colorbar; 
    axis xy;
    %imagesc(gridax,gridax,real(E(seti.meask2Ker(4,:)))); colorbar;
    title('changestructure')
end

if 0
    figure(104); %OK
    imagesc(real(Fq)); axis xy; colorbar; % colormap(litman); caxis([-1,3]); colorbar; % OK
    title('real(Fq)) - rebis changestructure')
end

if 0
    figure(105); % OK ?! (old remark: not when k^2 is not in kernel)
    gridax = seti.grid(2,1:seti.nCD);
    imagesc(gridax,gridax,real(E(DG))); axis xy; caxis([-0.8,+0.8]); colorbar; %colormap(litman); caxis([-1,0]); colorbar;
    title('DG - rebis changestructure')
end

if 0
    figure(106);
    imagesc(real(JA)); axis xy; colorbar; % colormap(litman); caxis([-1,3]); colorbar;
    title('real(JA)) - rebis changestructure')
    print('-f106','-dpng','figures/figure4');
    
    figure(107);
    imagesc(real(JB)); axis xy; colorbar; % colormap(litman); caxis([-1,3]); colorbar;
    title('real(JB) - rebis changestructure')
    print('-f107','-dpng','figures/figure6');
end

%%
% *Make test data --------------------------------------------------------*

if useSmoothTest == 0
    rng(0);
    h = randn(size(seti.qROIexact))+1i*randn(size(seti.qROIexact)); %test vector
else
    % -- use a smooth test vector (better convergence to zero...)
    mu1 = 0;
    sigma1 = 1/30;
    mu2 = seti.rCD/2;
    sigma2 = 1/40;
    %rROI = max(max(abs(seti.gridROI)));
    x = linspace(-seti.rCD,+seti.rCD,length(seti.qROIexact))'; %1...
    %plot(1:length(seti.qROIexact),pdf('Normal',x,mu,sigma));
    h = pdf('Normal',x,mu1,sigma1) + 1i*pdf('Normal',x,mu2,sigma2);
    % --
end

% h = reshape(h,seti.nROI,seti.nROI);

%%
% *Test derivative of FF(q) by Taylor-expansion ---------------------------*
DFFqh = JA*diagsparse(h)*JB; % FF'(q)[h]
if seti.dim == 2
    % Nice example to compare single and double precision in case of 
    % twoCornersOneBall2D (even nicer 1:8):
    tMat = 10.^-(1:6); 
else
    tMat = 10.^-(1:5);
end

dMat = size(tMat);
dMatSingle = size(tMat);
for i = 1:length(tMat)
    t = tMat(i);
    % FFqMeasph = [FF(q+t*h)] - 0
    FFqMeasph = mimo(seti, seti.qROIexact+t*h, 'simo');
    %fprintf('t = %e | max(max(FFqph)): %e | max(max(FFq)): %e | max(max(t*DFFqh)): %e \n',t,max(max(FFqph)),max(max(FFq)),max(max(t*DFFqh)));
    numerator = normws2(FFqMeasph - FFqMeas - t*DFFqh,seti);
    denominator = normroi2(t*h,seti);
    fprintf('t = %e | Numerator: %e | Denominator: %e | double prec. (default)\n',t,numerator, denominator);
    dMat(i) = normws2(FFqMeasph - FFqMeas - t*DFFqh,seti)/normroi2(t*h,seti);
    % Using that FF'(q)[.] is linear, so FF'(q)[t*h] = t*FF'(q)[h].
    % d = \|FF(q+t*h) - FF(q) - t*FF'(q)[h]\|_normws2 / \|t*h\|_normroi2

    if seti.dim == 2 % possible in 3D too but to save time only in 2D
        % -- test with single precision:
        FFqMeasph = mimo(seti, single(seti.qROIexact)+single(t)*single(h), 'simo');
        numerator = normws2(FFqMeasph - FFqMeas - t*DFFqh,seti);
        denominator = normroi2(t*h,seti);
        fprintf('t = %e | Numerator: %e | Denominator: %e | single prec.\n',t,numerator, denominator);
        dMatSingle(i) = normws2(FFqMeasph - FFqMeas - t*DFFqh,seti)/normroi2(t*h,seti);
        % --
    else
        dMatSingle(i) = 0;
    end
end

if seti.dim == 2
    figNo = 11;
elseif seti.dim == 3
    figNo = 15;
end
figure(figNo);
loglog(tMat,dMat,'ro',tMat,dMatSingle,'bo');
legend('double','single');
xlabel('t');
%ylabel('||FF(q+t*h) - FF(q) - FF''(q)[h]||_{WS} / ||t*h||_{ROI}');
ylabel('||FF(q+t*h) - FF(q) - FF''(q)[t*h]||_{WS} / ||t*h||_{ROI}');
%title('DG - rebis');
title('test derivative of FF(q)');
savePngFig(figNo,0,seti);

%%
% *Test Adjoint by inner products ----------------------------------------*

test = randn(size(seti.qROIexact))+1i*randn(size(seti.qROIexact));
fprintf('\n');
fprintf('Choose random vector g ...\n');
fprintf('In test case: y = 0 was choosen.\n');
s1 = innerroi(ADFFq,test,seti);
s2 = innerhs(FFqMeas,JA*diag(test)*JB,seti);
s = s1-s2;
fprintf('<FF''(q)^ast[FF(q)-y],g> - <FF(q)-y, FF''(q)[g]> = %g+i%g \n', ...
    real(s),imag(s));


%% 
% *Test adjoint for F(q) - y by direct computation -----------------------*
%
% See <ADFFqFast.html> for the method.

DG2 = zeros(size(ADFFq));
for iRX = 1:size(JA,1)
    for iTX = 1:size(JB,2)
        DG2 = DG2 + FFqMeas(iRX,iTX)*conj(JA(iRX,:).'.*JB(:,iTX));
    end;
end;
DG2 = diag(seti.dSMeas)/seti.dV*DG2;

% E = @(x) extendROItoCD(x,seti.ROImask);
% N = seti.nCD;

N = seti.nROI;
s = N*ones(1,seti.dim);
E = @(x) reshape(x, s);

EDG  = real(E(ADFFq));
EDG2 = real(E(DG2));

EDGi  = imag(E(ADFFq));
EDG2i = imag(E(DG2));

switch seti.dim
    case 2
        sliceWarn = '';
    case 3
        % Multiple slices by
        % v = round(linspace(1, N, min(10,N)), 0);
        % are not compatible with MATLAB -R2012a, so use:
        v = round(linspace(1, N, min(10,N)));
        
        figure(16); slice(EDG,v, N, 1); 
        title('Adjoint of derivative (ADFFq) computed by mimo'); 
        colormap(litman); colorbar;shading flat;
        
        figure(17); slice(EDG2,v, N, 1); 
        title('Adjoint of derivative (ADFFq) computed by direct A,B decomp.'); 
        colormap(litman); colorbar;shading flat;
        
        savePngFig(16,0,seti);
        savePngFig(17,0,seti);
        
        % Display a slice of the data:
        sliceWarn = ' (slice z = 0)';
        N2 = ceil((N-1)/2);
        EDG = EDG(:,:,N2);
        EDG2 = EDG2(:,:,N2);
end

if 0
    figure(108);imagesc(EDG);colorbar;
    axis xy;
    title(['computed by mimo' sliceWarn]);
    figure(109);imagesc(EDG2);colorbar;
    axis xy;
    title(['computed by direct A,B decomp.' sliceWarn]);

    figure(110);imagesc(EDGi);colorbar;
    axis xy;
    title(['computed by mimo (imag part)' sliceWarn]);
    figure(111);imagesc(EDG2i);colorbar;
    axis xy;
    title(['computed by direct A,B decomp. (imag part)' sliceWarn]);
end
