
if 0
    %% dense
    % parpool(4)    
    printFileName = 'solvePOPdense.csv';
    fid = fopen(printFileName, 'a');
    fprintf(fid,'degree,n,isBin,addComp,rOrder,solver,LBv,sec,iter,termcode\n');
    fclose(fid);
    
    clear paramset
    d = 3;
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {true}; %%{false, true}; %comp
    paramset{3} = {6}; %% {30};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    
    params = genAllComb(paramset);
    
    binset  = params(:, 1);
    compset = params(:, 2);
    sizeset = params(:, 3);
    solset  = params(:, 4);
    
    % parfor ii = 1:length(binset)
    for ii = 1:length(binset)
        % solvePOPdense(d,n,isBin,addComplement,solver)
        % solvePOPdense(d, sizeset{ii}, binset{ii}, compset{ii},solset{ii});
        solvePOPdense(printFileName,d, sizeset{ii}, binset{ii}, compset{ii}, solset{ii});
    end
end

if 0
    %% arrow
    % parpool(4)    
    printFileName = 'solvePOParrow.csv';
    fid = fopen(printFileName, 'a');
    fprintf(fid,'degree,n,a,b,c,l,isBin,addComp,rOrder,solver,LBv,sec,iter(APGR),termcode,iterBP\n');
    close(fid);
    
    clear paramset
    d = 5; a = 10;
    paramset{1} = {true}; %% {false}; %% {true, false}; %bin box
    paramset{2} = {true}; %% {false, true}; %comp
    paramset{3} = {2}; %% {2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    
    params = genAllComb(paramset);
    
    binset  = params(:,1);
    compset = params(:,2);
    sizeset = params(:,3);
    solset  = params(:,4);
    
    % parfor ii = 1:length(binset)
    for ii = 1:length(binset)
        %solvePOParrow(d,a,b,c,l,isBin,addComplement);
        solvePOParrow(printFileName,d,a,2,2,sizeset{ii},binset{ii},compset{ii},solset{ii});
    end
end

if 0
    %% chordal
    % parpool(4)    
    printFileName = 'solvePOPchordal.csv';
    fid = fopen(printFileName, 'a');
    fprintf(fid,'degree,n,rRange,isBin,addComp,rOder,solver,LBv,sec,iter(APGR),termcode,iterBP\n');    
    close(fid);
    
    clear paramset
    d = 3; %% 5
    paramset{1} = {true}; % bin box
    paramset{2} = {false}; %% {false, true}; %comp
    paramset{3} = {100}; %% {240}; %size
    paramset{4} = {'BP'}; %% {'BP','sdpnal'};
    
    params = genAllComb(paramset);
    
    binset = params(:,1);
    compset = params(:,2);
    sizeset = params(:,3);
    solset = params(:,4);
    
    % parfor ii = 1:length(binset)
    for ii = 1:length(binset)
        radiorange = 0.1; 
        solvePOPchordal(printFileName,d,sizeset{ii},radiorange,binset{ii},compset{ii},solset{ii});
    end
end

if 0
    %% dense
    % parpool(4)    
    printFileName = 'solvePOPdense.csv';
    fid = fopen(printFileName, 'a');
    fprintf(fid,'degree,n,isBin,addComp,rOrder,solver,LBv,sec,iter,termcode\n');
    fclose(fid);
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {800,1200,1600};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {2}; %d
    
    params2 = genAllComb(paramset);
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {40,60};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {3}; %d
    
    params3 = genAllComb(paramset);
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {10,15};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {5}; %d
    
    params5 = genAllComb(paramset);
    
    binset  = [params2(:, 1); params3(:, 1); params5(:, 1)];
    compset = [params2(:, 2); params3(:, 2); params5(:, 2)];
    sizeset = [params2(:, 3); params3(:, 3); params5(:, 3)];
    solset  = [params2(:, 4); params3(:, 4); params5(:, 4)];
    dset    = [params2(:, 5); params3(:, 5); params5(:, 5)];
    
    parfor ii = 1:length(binset)
        solvePOPdense(printFileName,dset{ii}, sizeset{ii}, binset{ii}, compset{ii},solset{ii});
    end
end

if 0
    %% arrow
    % parpool(4)    
    printFileName = 'solvePOParrow.csv';
    fid = fopen(printFileName, 'a');
    fprintf(fid,'degree,n,a,b,c,l,isBin,addComp,rOrder,solver,LBv,sec,iter(APGR),termcode,iterBP\n');
    close(fid);

    a = 5;
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {40,60};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {6}; %d
    
    params6t = genAllComb(paramset);
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {40,60};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {8}; %d
    
    params8t = genAllComb(paramset);
    
    clear paramset
    paramset{1} = {false}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {16,20};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {6}; %d
    
    params6f = genAllComb(paramset);
    
    clear paramset
    paramset{1} = {false}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {9};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {8}; %d
    
    params8f = genAllComb(paramset);
    
    
    binset  = [params6t(:, 1); params8t(:, 1); params6f(:, 1); params8f(:, 1)];
    compset = [params6t(:, 2); params8t(:, 2); params6f(:, 2); params8f(:, 2)];
    sizeset = [params6t(:, 3); params8t(:, 3); params6f(:, 3); params8f(:, 3)];
    solset  = [params6t(:, 4); params8t(:, 4); params6f(:, 4); params8f(:, 4)];
    dset    = [params6t(:, 5); params8t(:, 5); params6f(:, 5); params8f(:, 5)];
    
    parfor ii = 1:length(binset)
        solvePOParrow(printFileName,dset{ii},a,2,2,sizeset{ii},binset{ii},compset{ii},solset{ii});
    end
end

if 0
    %% dense
    % parpool(4)
    printFileName = 'solvePOPdense.csv';
    fid = fopen(printFileName, 'a');
    fprintf(fid,'degree,n,isBin,addComp,rOrder,solver,LBv,sec,iter,termcode\n');
    fclose(fid);
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {50};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {3}; %d
    
    params3 = genAllComb(paramset);
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {20};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {5}; %d
    
    params5 = genAllComb(paramset);
    
    binset  = [ params3(:, 1); params5(:, 1)];
    compset = [ params3(:, 2); params5(:, 2)];
    sizeset = [ params3(:, 3); params5(:, 3)];
    solset  = [ params3(:, 4); params5(:, 4)];
    dset    = [ params3(:, 5); params5(:, 5)];
    
    parfor ii = 1:length(binset)
        solvePOPdense(printFileName,dset{ii}, sizeset{ii}, binset{ii}, compset{ii}, solset{ii});
    end
end

if 0
    %% arrow
    % parpool(4)
    printFileName = 'solvePOParrow.csv';
    fid = fopen(printFileName, 'a');
    fprintf(fid,'degree,n,a,b,c,l,isBin,addComp,rOrder,solver,LBv,sec,iter(APGR),termcode,iterBP\n');
    close(fid);
    
    a = 5;
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {80, 100};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {6}; %d
    
    params6t = genAllComb(paramset);
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {80};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {8}; %d
    
    params8t = genAllComb(paramset);
    
    binset  = [params6t(:, 1); params8t(:, 1)];
    compset = [params6t(:, 2); params8t(:, 2)];
    sizeset = [params6t(:, 3); params8t(:, 3)];
    solset  = [params6t(:, 4); params8t(:, 4)];
    dset    = [params6t(:, 5); params8t(:, 5)];
    
    parfor ii = 1:length(binset)
        solvePOParrow(printFileName,dset{ii},a,2,2,sizeset{ii},binset{ii},compset{ii},solset{ii});
    end
end

if 0
    %% chordal
    % parpool(4)
    printFileName = 'solvePOPchordal.csv';
    fid = fopen(printFileName, 'a');
    fprintf(fid,'degree,n,rRange,isBin,addComp,rOder,solver,LBv,sec,iter(APGR),termcode,iterBP\n');    
    close(fid);
    
    clear paramset
    paramset{1} = {false}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {140};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {6,8}; %d
    
    params1 = genAllComb(paramset);
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {170,230};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {6}; %d
    
    params2 = genAllComb(paramset);
    
    clear paramset
    paramset{1} = {true}; %%, false}; %bin box
    paramset{2} = {false, true}; %comp
    paramset{3} = {190};%{2, 3, 4}; %size
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    paramset{5} = {8}; %d
    
    params3 = genAllComb(paramset);
    
    binset  = [params1(:, 1); params2(:, 1); params3(:, 1)];
    compset = [params1(:, 2); params2(:, 2); params3(:, 2)];
    paramset{4} = {'BP'}; %% {'sdpnal', 'BP'};
    dset    = [params1(:, 5); params2(:, 5); params3(:, 5)];
    
    parfor ii = 1:length(binset)
        solvePOPchordal(printFileName,dset{ii},sizeset{ii},0.1,binset{ii},compset{ii},solset{ii});
    end
end

