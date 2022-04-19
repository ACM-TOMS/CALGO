init;
seti.expData = 'fresnel';
seti.fresnelFreq = 5*1E9;
seti.fresnelFile = 'inexpdata/fresnel_opus_1/twodielTM_8f.exp';
seti = setData(seti);
imagesc(real(seti.G(seti.incField(:,6)))); colorbar; axis xy;
