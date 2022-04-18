% Example to the guide: Working with data from Institute Fresnel

init;
filename = 'inexpdata/fresnel_opus_1/twodielTM_8f.exp';
[uTotRX, uIncRX, frequencies, rTX, nTX, rRX, nRX] = readRAWData(filename);
% figure(1); imagesc(real(uIncRX(:,:,1))); colorbar; axis xy;

% Publication: plot such that NaN is easily visible and save
colormap default; cmap = colormap; cmap(1,:) = 1;
data = real(uIncRX(:,:,1));
data(isnan(data)) = min(min(data))-0.1;
figure(1); imagesc(data); colormap(cmap); colorbar; axis xy;
set(gca,'FontSize',20); axis square; print(1,'-depsc','guideWorkingFresnel.eps');
