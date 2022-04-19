init;
[seti,uTotRX,uIncRX,uScaRX] = matchIncFieldTrans('inexpdata/fresnel_opus_1/twodielTM_8f.exp',5*1E9);

% Compute the incident field uIncROI on ROI
% seti.ampCalc = 1;   % method to compute the coefficients
seti.nuMax = 7;     % 2*nuMax+1 coefficients are computed for polynomial approximation
[uIncROI,errC] = matchIncField(uIncRX,seti,'ROI');

figure(1); imagesc(real(reshape(uIncROI(:,6),[seti.nROI seti.nROI]))); colorbar; axis xy;
% Note that uIncROI(:,n) reads the incident field of the n-th transmitter

%% Save plot
%set(gca,'FontSize',20); axis square; print(1,'-depsc','guideMatchIncField.eps');
