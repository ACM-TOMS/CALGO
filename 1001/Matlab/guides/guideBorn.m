init;

% k: wave number
% k = 5;
k = 250;

%%
seti.k = k;

seti.contrast = 'corner2D';
seti = setData(seti);
qROI = seti.qROIexact;

uIncROI = seti.dSInc.*seti.incField(:,1); % incident field of first transmitter (1, ..., seti.incNb)

%% uBorn

% Solution 1:
% Define factor k^2 in the volume potential operator V.
V  = @(x) seti.k^2.*helmholtz2Dr2r(x, seti); % volume potential operator V: L^2(ROI) -> L^2(ROI)
QU = @(x) qROI .* x;
uBorn1 = V(QU(uIncROI));

% Solution 2:
% As solution 1, but V and QU are defined in intOpsFuncs.
% [V, ~, ~, ~, ~, QU, ~, ~, ~] = intOpsFuncs(qROI, seti);
% uBorn2 = V(QU(uIncROI));

% Solution 3:
% uBorn3 = seti.k^2.*helmholtz2Dr2r(qROI.*uIncROI, seti);

% Test:
% norm(uBorn1-uBorn2)
% norm(uBorn2-uBorn3)

%% uScattROI

% Solution 1a:
% uScattROI1 = solveLippmannSchwinger(@(x) V(QU(x)), V(QU(uIncROI)), seti); 

% Solution 1b:
uScattROI1 = solveLippmannSchwinger(@(x) V(qROI.*x), V(qROI.*uIncROI), seti); 

% Solution 2:
% [~, uScattROI2] = simo(qROI,uIncROI,seti);

% Test:
% norm(uScattROI1-uScattROI2)

%% uScattROI and uBorn are similar in case of small wavelengths

rel = norm(uBorn1-uScattROI1)/norm(uScattROI1);
fprintf('  Wave number: %g\n',seti.k);
fprintf('  Relative error: %g\n',rel');

%%
% 
% figure(1); f1 = imagesc(real(seti.G(qROI))); colorbar; axis xy; % predefined contrast
% set(gca,'FontSize',20); axis square; print(1,'-depsc','pi3born1.eps');
% figure(2); f2 = imagesc(real(seti.G(uIncROI))); colorbar; axis xy; % incident field
% set(gca,'FontSize',20); axis square; print(2,'-depsc','pi3born2.eps');
% figure(3); f3 = imagesc(real(seti.G(uBorn1))); colorbar; axis xy; % Born approximation of scattered field
% set(gca,'FontSize',20); axis square; print(3,'-depsc','pi3born3.eps');
% figure(4); f4 = imagesc(real(seti.G(uScattROI1))); colorbar; axis xy; % Scattered field (by solving Lippmann Schwinger eq.)
% set(gca,'FontSize',20); axis square; print(4,'-depsc','pi3born4.eps');
