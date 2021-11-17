% runme.m
%{
	*** MATLAB script file to build and run a sample program using
	***            Talbot Suite and Talbot Suite DE
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
OUTdir = './';                        % OUTPUT directory
TALBOTdir = '../../../TalbotSuiteDE'; % Talbot Suite DE directory
EXdir = '../../';                     % Example 0b (Duffy) directory

%% BUILD EXECUTABLE
% shared functions
FILENAMES = ['./num_sol.c ./MEX_complexArray.c ' EXdir '/LTsings2.c ' TALBOTdir '/COM/COM_Talbot_pack.c ' TALBOTdir '/COM_DE/COM_Talbot_pack_DE.c '];

% for Talbot Suite (copied in current folder with no correction to NOPTS)
FILENAMES = [FILENAMES EXdir '/LTfun.c ./SEQ_Talbot_pack.c '];

% for Talbot Suite DE
FILENAMES = [FILENAMES './SEQ_LTsamples_fun.c ' EXdir '/LTfun2.c ' TALBOTdir '/FUN_DE/SEQ_Talbot_pack_DE.c '];

fprintf('\nBUILD\n\n')
% TO COMPILE FOR ACCURACY (option -v: verbose)
STRING = ['mex ' OPTS FILENAMES];
%STRING = ['mex ' OPTS ['./num_sol.c ' EXdir '/LTfun.c ' EXdir '/ILTfun2.c ' FILENAMES]];

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
