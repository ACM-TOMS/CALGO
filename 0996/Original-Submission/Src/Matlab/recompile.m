dirlist = cell(0);
filelist = cell(0);
dirlist{end+1} = 'relaxation'; filelist{end+1} = 'colgcdsp.cpp';
dirlist{end+1} = 'relaxation'; filelist{end+1} = 'cutByDinic.cpp';
dirlist{end+1} = 'polyTools'; filelist{end+1} = 'enumCombWithNoRep.cpp';
dirlist{end+1} = 'polyTools'; filelist{end+1} = 'enumCombWithRep.cpp';
dirlist{end+1} = 'solver/BP'; filelist{end+1} = 'mexProxMonotonicNew.cpp';
dirlist{end+1} = 'solver/BP/mexfun_TohK'; filelist{end+1} = 'blkEigMex.c';
dirlist{end+1} = 'solver/BP/mexfun_TohK'; filelist{end+1} = 'mexeigK.c';
dirlist{end+1} = 'solver/BP/mexfun_TohK'; filelist{end+1} = 'mexFnormK.c';
dirlist{end+1} = 'solver/BP/mexfun_TohK'; filelist{end+1} = 'mexeigPartialK.c';

fprintf('List: ');
for ii = 1:length(filelist)
    fprintf(' %s', filelist{ii});
end
fprintf('\n');

for ii = 1:length(filelist)
    fprintf(['Compiling ' dirlist{ii} '/' filelist{ii} '...\n'])
    eval(['mex -O -largeArrayDims ' dirlist{ii} '/' filelist{ii} ' -outdir ' dirlist{ii} ' -lmwlapack -lmwblas']);
end
fprintf('DNNPOP is successfully installed!\n');
