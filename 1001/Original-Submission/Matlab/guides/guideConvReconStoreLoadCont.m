seti.tau = 2.0;
seti.loadqROIcomp = sprintf('%s/save_qROIcomp_iOutStop.mat',seti.dirname);
seti = setRecon(seti,4,2);
seti = recon(seti,4,2);
