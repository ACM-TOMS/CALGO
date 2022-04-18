% runme.m
%{
	*** MATLAB script file to build and run a sample program using
	***              Talbot Suite DE
	***
	*** author: Mariarosaria Rizzardi
	***         mariarosaria.rizzardi@uniparthenope.it
	***
%}
clear; clc

%% Check for the platform
%{
% computer = PCWIN              MATLAB and Windows 32 bit
%            PCWIN64            MATLAB and Windows 64 bit
%            GLNX86             MATLAB and Unix 32 bit
%            GLNXA64            MATLAB and Unix 64 bit
disp(['MATLAB is running on this platform: ' computer])
if ispc       % MATLAB + Windows
    OPTS = ' COMPFLAGS="$COMPFLAGS -fopenmp" LINKFLAGS="$LINKFLAGS -fopenmp -pthread" LINK_LIB="$LINK_LIB -lgomp" '; % with OpenMP
elseif isunix % MATLAB + Unix
    OPTS = ' ''CFLAGS="\$CFLAGS -fopenmp"'' -lgomp -lrt ';              % Link with time library
elseif ismac % MATLAB + MAC
    error('For MAC platform check options');
else
    error('Cannot check the running platform');
end
%}
OPTS = [];

%% DIRECTORIES
OUTdir = './';                         % OUTPUT directory
TALBOTdir = '../../../TalbotSuiteDE';  % Talbot Suite DE directory
EXdir = '../..';                       % Example 3a directory

%% SHARED FILES
FILENAMES = './SEQ_LTsamples_ode45.c ./MEX_complexArray.c';
FILENAMES = [FILENAMES ' ' EXdir '/LTsings2.c'];
FILENAMES = [FILENAMES ' ' TALBOTdir '/COM/COM_Talbot_pack.c ' TALBOTdir '/COM_DE/COM_Talbot_pack_DE.c ' TALBOTdir '/FUN_DE/SEQ_Talbot_pack_DE.c'];


%% BUILD EXECUTABLE
choice = menu('Select:','1)  accuracy test','2)  time test');

switch choice
    case 1
        fprintf('entered 1: accuracy test\n')
        fprintf('\nBUILD\n\n')
        % TO COMPILE FOR ACCURACY (option -v: verbose)
        STRING = ['mex ' OPTS ['./rel_err.c ' EXdir '/ILTfun2.c ' FILENAMES]];
    otherwise
        fprintf('entered 2: time test\n')
        fprintf('\nBUILD\n\n')
        % TO COMPILE FOR TIMES (option -v: verbose)
        STRING = ['mex ' OPTS ['./SEQ_PDE_skill_TIMES.c ' FILENAMES]];
end

fprintf('\nExecute:\n'); disp(['>> ' STRING]);
eval(STRING)
fprintf('Done.\n')
fprintf('\n********************************************\n\n')

%% RUN EXECUTABLE
fprintf('RUN EXECUTABLE\n\tMAIN\n')
MAIN
fprintf('\n\n********************************************\n')
fprintf('\tTO REMOVE EXECUTABLE FILES,\n')
fprintf('\tENTER:\n')
fprintf('\t\tclean\n\n')
fprintf('\tOTHERWISE RUN AS:\n')
fprintf('\t\t[clear]  optional\n')
fprintf('\t\tMAIN')
fprintf('\n\n********************************************\n\n')
