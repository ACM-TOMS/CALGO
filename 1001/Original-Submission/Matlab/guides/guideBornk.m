N = 100; % N = 70; % N = 20;
k = zeros(1,N);
rel = zeros(1,N);
for i = 1:N
  k(i) = i*1;
  seti.k = k(i);
  
  seti = setData(seti);
  qROI = seti.qROIexact;
  uIncROI = seti.dSInc.*seti.incField(:,1); % incident field of first transmitter (1, ..., seti.incNb)

  V  = @(x) seti.k^2.*helmholtz2Dr2r(x, seti); % volume potential operator V: L^2(ROI) -> L^2(ROI)
  QU = @(x) qROI .* x;
  uBorn1 = V(QU(uIncROI));

  uScattROI1 = solveLippmannSchwinger(@(x) V(QU(x)), V(QU(uIncROI)), seti); 

  % Relative error of Born approximation uBorn compared with scattered field uScattROI
  rel(i) = norm(uBorn1-uScattROI1)/norm(uScattROI1);

end

% uScattROI and uBorn are similar in case of small wavelengths
% Plot the relative error "rel" in dependence of the wwave number "k"
figure(1); h = plot(k,rel);
set(gca,'FontSize',20); set(h,'LineWidth',2); axis square; print(1,'-depsc','guideBornk.eps');
