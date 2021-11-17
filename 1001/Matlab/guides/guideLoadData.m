init;
seti.expData = 'fresnel';
seti.rCD = 0.2;
seti.nCD = 256;
seti.fresnelFreq = 5*1E9;
seti.fresnelFile = 'inexpdata/fresnel_opus_1/twodielTM_8f.exp';
% seti.ampCalc = 1;   % method to compute the coefficients
seti.nuMax = 7;     % 2*nuMax+1 coefficients are computed for polynomial approximation
seti = checkConsisExpData(seti);
seti = loadData(seti);

%%

figure(1); imagesc(real(seti.G(seti.incField(:,6)))); colorbar; axis xy;
% Note that seti.incField(:,n) reads the incident field of the n-th transmitter.

%% Save plot
%set(gca,'FontSize',20); axis square; print(1,'-depsc','guideLoadData.eps');
