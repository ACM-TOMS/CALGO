init;
for N = [1000] % N = [50, 100, 500, 1000, 5000, 10000]
    ticBorn = tic;
    seti.nOut = 1;
    seti.pdaN = N;
    seti.k = 70;
    disp(' ');
    fprintf('pdaN = %i.\n',seti.pdaN);
    seti = setData(seti);
    seti = setRecon(seti);
    seti = recon(seti);
    
    tocBorn = toc(ticBorn);
    fprintf('Elapsed time is %05.1f min.\n',tocBorn/60);

    figure(1); imagesc(real(seti.G(seti.qROIcomp))); colorbar; axis xy;
    set(gca,'FontSize',20); axis square; caxis([-0.15 +1.42]); 
    filename = sprintf('guideBornInv_%i.eps',seti.pdaN);
    print(1,'-depsc',filename);

    fprintf('rel. dis. %.4f...\n', seti.dis(seti.iOutStop));
    fprintf('rel. err. %.4f.\n', seti.err(seti.iOutStop));
end

