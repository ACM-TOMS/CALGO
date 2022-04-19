% build_accuracy.m
% Compile, with MATLAB mex, and build mex-executable
%	***   MATLAB + Unix/Windows version   ***

%% Check for the platform [OpenMP option]
%{
% computer = PCWIN              MATLAB and Windows 32 bit
%            PCWIN64            MATLAB and Windows 64 bit
%            GLNX86             MATLAB and Unix 32 bit
%            GLNXA64            MATLAB and Unix 64 bit
disp(['MATLAB is running on this platform: ' computer])
if ispc       % MATLAB + Windows
    OPTS = ' COMPFLAGS="$COMPFLAGS -fopenmp" LINKFLAGS="$LINKFLAGS -fopenmp -pthread" LINK_LIB="$LINK_LIB -lgomp" '; % with OpenMP
elseif isunix % MATLAB + Unix
    OPTS = ' ''CFLAGS="\$CFLAGS -fopenmp"'' -lgomp '; % Link without time library
elseif ismac % MATLAB + MAC
    error('For MAC platform check options');
else
    error('Cannot check the running platform');
end
%}
OPTS = [];

%% DIRECTORIES
OUTdir = './';                        % OUTPUT directory
TALBOTdir = '../../../TalbotSuiteDE'; % Talbot Suite DE directory
EXdir  = '../..';                     % Example 5 directory

%% BUILD EXECUTABLE
FILENAMES = './SEQ_main_ACCURACY.c ./SEQ_LTsamples_blktrd.c ./MEX_complexArray.c';
FILENAMES = [FILENAMES ' ' EXdir '/ILTfun2.c ' EXdir '/LTsings2.c'];
FILENAMES = [FILENAMES ' ' TALBOTdir '/COM/COM_Talbot_pack.c ' TALBOTdir '/COM_DE/COM_Talbot_pack_DE.c ' TALBOTdir '/FUN_DE/SEQ_Talbot_pack_DE.c'];

%% TO COMPILE (option -v: verbose)
STRING = ['mex ' OPTS FILENAMES];
fprintf('\nExecute:\n'); disp(['>> ' STRING]);
eval(STRING)
fprintf('Done.\n')

%% TO RUN: SEQ_main_ACCURACY for default tolerance or ...
fprintf('\nTO RUN EXECUTE:')
fprintf('\n\t>> tol=1e-6; SEQ_main_ACCURACY(tol)\n')
