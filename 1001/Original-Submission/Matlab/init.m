%% init
% Initialization: add required pathes.
%
%% Syntax
%
%   init;
%
%% Description
%
% * Add required pathes for functions
% * Optional: set parameter "closed" (0 or 1)
%   (1 does mean: closed code is available). 
%   In public version |closed = 0|. 
%   (If closed is not set it is set to 0 in |setInput| automatically.)
%
%% More About
%
% New developments are stored in the directory "closed" and other 
% directoris with suffix "Closed". These folders are not published in this 
% version of code. Therefore we introduce the variable
% "closed".
% In published code it will be 0, because the directory does not exist.
%
% The value of |closed| will be stored in |seti.closed| in |setInput.m|
% (afterwards |closed| is cleared).
%
%% See Also
%
% * <start.html>
%
%% Code

%%
% *input*

addpath incontrasts/           % input contrasts
addpath incontrasts/2D/        % 2D
addpath incontrasts/2DFresnel/ % 2D compatible to real-world data from Institute Fresnel
addpath incontrasts/3D/        % 3D

% addpath inseti/         % input parameter settings
addpath(genpath('inseti/')) % add subfolders of inseti too

%%
% *process*

addpath proc/           % process in general
addpath proc/auxi/      % auxiliary functions (Do not rename this folder name in "aux", because Microsoft Windows can not deal with such a name! (Linux is no problem.))
addpath proc/expData/   % functions to deal with real-world (experimentally measured) data from Institute Fresnel
addpath proc/expSetup/  % functions to deal with experimental set-up (positions of transmitters and receivers; incident )
addpath proc/intOps/    % intergral operators
addpath proc/norms/     % definitions of norms
addpath proc/operators/ % operators
addpath proc/plots/     % plots
addpath proc/plotsAux/  % auxilary functions to plot (save, colormap, ...)
addpath proc/setData/   % setting of geometry and simulation
addpath proc/setInput/  % process of general input (make directories...)
addpath proc/recon/     % reconstruction process
addpath proc/reconAux/  % auxiliary files for reconstruction

%%
% *3rdparty*

addpath 3rdparty/ 	% 3rd party code

%% 
% *tests*

addpath tests/          % test functions
addpath tests/auxi/      % auxiliary functions (e.g. reference data)

%%
% *convenience functions*
%
% Some convenience functions (are not used in other parts of program...)

addpath conv/

%%
% *guides*
%
% Some convenience functions for the guides (are not used in other parts of program...)

addpath guides/
addpath guides/auxi/

%%
% *data structure reference*
%
% Add the path such that content of folder "dataStructure" is considered
% when docCreate is used.

addpath docCreate/addMfiles/

%%
% *Option: closed code*
%
% This part can be commented out (or set closed = 0).
% 
closed = 0; % 0 or 1 (1 does mean: closed code is available)

if closed == 1
disp(' ')
disp('################ --- Closed code is available ---- ################')
disp(' ');
end

if closed == 1

    % closed
    addpath closed/
    
    % input
    addpath incontrasts/2DClosed/
    addpath incontrasts/3DClosed/
    % addpath insetiClosed/
    addpath(genpath('insetiClosed/')) % add subfolders of insetiClosed too
    
    % process
    addpath proc/auxiClosed/
    addpath proc/expDataClosed/
    addpath proc/intOpsClosed/
    addpath proc/operatorsClosed/
    addpath proc/plotsClosed/
    addpath proc/reconAuxClosed/
    addpath proc/reconClosed/
    addpath proc/setDataClosed/
    
    % tests
    addpath testsClosed/
end
