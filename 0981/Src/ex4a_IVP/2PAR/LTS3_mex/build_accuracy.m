% build_accuracy.m
% Compile, with MATLAB mex, and build mex-executable
%	***   MATLAB + Unix/Windows version   ***

%% CHECK for Parallel Computing Toolbox
fun = 'parpool';      % ==>   distcomp
pat = '(?<=^.+[\\/]toolbox[\\/])[^\\/]+'; % match regular expression (case sensitive)
if isempty(regexp(which(fun), pat, 'match', 'once'))
    error('The Parallel Computing Toolbox is required to run this program')
else
    fprintf('\n*** The Parallel Computing Toolbox is installed. ***\n')
end

%% Check for the platform
% {
% computer = PCWIN              MATLAB and Windows 32 bit
%            PCWIN64            MATLAB and Windows 64 bit
%            GLNX86             MATLAB and Unix 32 bit
%            GLNXA64            MATLAB and Unix 64 bit
disp(['MATLAB is running on this platform: ' computer])
if ispc       % MATLAB + Windows
    OPTS = ' COMPFLAGS="$COMPFLAGS -fopenmp" LINKFLAGS="$LINKFLAGS -fopenmp -pthread" LINK_LIB="$LINK_LIB -lgomp" '; % with OpenMP
elseif isunix % MATLAB + Unix
    OPTS = ' ''CFLAGS="\$CFLAGS -fopenmp"'' -lgomp ';                   % Link without time library
elseif ismac % MATLAB + MAC
    error('For MAC platform check options');
else
    error('Cannot check the running platform');
end
%}
%OPTS = [];

%% DIRECTORIES
OUTdir = './';                        % OUTPUT directory
TALBOTdir = '../../../TalbotSuiteDE'; % Talbot Suite DE directory
EXdir  = '../..';                     % Example 4a directory

%% BUILD EXECUTABLE
FILENAMES = './OMP_main_ACCURACY.c ./OMP_LTsamples_ode45.c ./OMP_MEX_complexArray.c';
FILENAMES = [FILENAMES ' ' EXdir '/ILTfun2.c ' EXdir '/LTsings2.c'];
FILENAMES = [FILENAMES ' ' TALBOTdir '/COM/COM_Talbot_pack.c ' TALBOTdir '/COM_DE/COM_Talbot_pack_DE.c ' TALBOTdir '/FUN_DE/OMP_Talbot_pack_DE.c'];

%% TO COMPILE (option -v: verbose)
STRING = ['mex ' OPTS FILENAMES];
fprintf('\nExecute:\n'); disp(['>> ' STRING]);
eval(STRING)
fprintf('Done.\n')

%% APPEND the "MinGW-w64\...\bin" folder to your system PATH variable in Windows
append_MinGW_dir

%% TO RUN: OMP_main_ACCURACY for default tolerance or ...
fprintf('\nTO RUN EXECUTE:')
fprintf('\n\t>> tol=1e-12; thrds1=4; OMP_main_ACCURACY(tol,thrds1)')
fprintf('\nOR')
fprintf('\n\t>> tol=1e-12; thrds1=4; thrds2=3; OMP_main_ACCURACY(tol,thrds1,thrds2)\n')

