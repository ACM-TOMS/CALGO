% runme.m
% Compile, with MATLAB mex, and build mex-executable
%	***   MATLAB + Unix/Windows version   ***

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
    OPTS = ' CLIBS="\$CLIBS -lrt" '; % Link with time library
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

%% BUILD EXECUTABLES
FILENAMES = ['./SEQ_Talbot1_SUM_mex.c ' TALBOTdir '/COM/COM_Talbot_pack.c ' TALBOTdir '/COM_DE/COM_Talbot_pack_DE.c ' TALBOTdir '/FUN_DE/SEQ_Talbot_pack_DE.c'];
STRING = ['mex ' OPTS FILENAMES]; % (option -v: verbose)
fprintf('\nExecute:\n'); disp(['>> ' STRING]);
eval(STRING)
fprintf('Done.\n')

%% TO RUN:
choice = menu('Select:','1)  accuracy test','2)  time test');

switch choice
    case 1
        fprintf('\nFor accuracy test, execute in MATLAB:\n\t>> MAIN_accuracy\n')
        MAIN_accuracy

    case 2
        fprintf('\nFor  time test,   execute in MATLAB:\n\t>> MAIN_times\n')
        MAIN_times
end



