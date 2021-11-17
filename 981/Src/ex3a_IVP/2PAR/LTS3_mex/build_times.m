% build_times.m
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
    %OPTS = ' CFLAGS="\$CFLAGS -fopenmp" CLIBS="\$CLIBS -lgomp -lrt" '; % Link with time library
    %OPTS = ' ''CFLAGS="\$CFLAGS -fopenmp"'' -lgomp ';                  % Link without time library
    OPTS = ' ''CFLAGS="\$CFLAGS -fopenmp"'' -lgomp -lrt ';              % Link with time library
elseif ismac % MATLAB + MAC
    error('For MAC platform check options');
    %OPTS = 'CFLAGS="\$CFLAGS -fopenmp" CLIBS="\$CLIBS -lgomp" '; % Link without time library
    %OPTS = 'CFLAGS="\$CFLAGS -fopenmp" CLIBS="\$CLIBS -lgomp -lrt" '; % Link with time library
else
    error('Cannot check the running platform');
end
%}
%OPTS=[];

%% DIRECTORIES
OUTdir = './';                        % OUTPUT directory
TALBOTdir = '../../../TalbotSuiteDE'; % Talbot Suite DE directory
EXdir  = '../..';                     % Example 3a directory

%% BUILD EXECUTABLE
FILENAMES = './OMP_main_TIMES.c ./OMP_LTsamples_ode45.c ./OMP_MEX_complexArray.c';
FILENAMES = [FILENAMES ' ' EXdir '/LTsings2.c'];
FILENAMES = [FILENAMES ' ' TALBOTdir '/COM/COM_Talbot_pack.c ' TALBOTdir '/COM_DE/COM_Talbot_pack_DE.c ' TALBOTdir '/FUN_DE/OMP_Talbot_pack_DE.c'];

%% TO COMPILE (option -v: verbose)
STRING = ['mex ' OPTS FILENAMES];
fprintf('\nExecute:\n'); disp(['>> ' STRING]);
eval(STRING)
fprintf('Done.\n')

%% APPEND the "MinGW-w64\...\bin" folder to your system PATH variable in Windows
append_MinGW_dir

%% TO RUN: OMP_main_TIMES for default tolerance or ...
fprintf('\nTO RUN EXECUTE:')
fprintf('\n\t>> tol=1e-6; jFUN=1; NTval=20; NXval=20; OMP_main_TIMES(tol,jFUN,NTval,NXval)\n')

