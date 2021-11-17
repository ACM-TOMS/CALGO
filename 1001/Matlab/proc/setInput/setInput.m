%% setInput
% Set input parameters.
%
% Warning: This function clears variables and creates a directory. Also the
% routine |eval| is used. Be careful when using.
%
%% Syntax
%
%   setInput;
%
%% Description
%
% |setInput| is called in <start.html>.
%
% The code in this file is separated in the following parts:
%
% * Various 1
% * Set parameters
% * Various 2
% * Set and make directories for output
% * File suffix
%
% We will give descriptions and examples for each part.
%
% If the mentioned arguments does not exist, they are set automatically in setInput.
%
% The central variable is the structural array |seti|.
%

%% Various 1

%%
% *Description of Various 1: deal with various input parameters*
%
% Use |varalpha|, |varbeta|, |vardelta| to call |start| with
% several inputs for regularization parameters $\alpha$, $\beta$, or noiselevels $\delta$.
%
% A parameter is set in the correspong file as in the following list:
%
%   usevaralpha = 1; % is set in varalpha
%   usevarbeta  = 1; % is set in varbeta
%   usedelta    = 1; % is set in vardelta
%
% If they are not set, they will be set to |0|.

%%
% *Example for Various 1: Example to deal with various input*
%
% To tackle the inverse scattering problem with several regularization
% parameters $\alpha$ do the following.
%
% Define the alpha-values in |varalpha.m|, e.g.:
%
%   alpha = [8E2; 2E3; 3E3];
%
% Use MATLAB commands to set parameters and contrast and call |varalpha|:
%
%   inseti = 'example';
%   varalpha;
%
% Results are stored in a new folder in |output|.

%%
% *Input Arguments of Various 1:*
%
% Note that all inputs are optional.
%
% * usevaralpha     :   If |varalpha| calls |start|, it is set to 1 (otherwise 0)
%                       and various alpha-values are used
%                       (alpha: regularization parameter of the sparsity penalty).
% * usevarbeta      :   If |varbeta| calls |start|, it is set to 1 (otherwise 0) 
%                       and various beta-values are used
%                       (beta: regularization parameter of the total variation penalty).
% * usevardelta     :   If |vardelta| calls |start|, it is set to 1 (otherwise 0) 
%                       and various relative noise levels delta are used.
%
%%
% *Output Arguments of Various 1*
%
% If mentioned input parameters are not set, they are defined and set to 0.
% Please be aware of
%
%   clearvars -except closed test inseti
%
% * seti.dirDatetime    : date and time (format: YYYYMMDDThhmmss)
% * seti.fileSuffix     : suffix name for filename
%    
%%
% *Code of Various 1*

if ~exist('usevaralpha','var')
    usevaralpha = 0;
end
if ~exist('usevarbeta','var')
    usevarbeta = 0;
end
if ~exist('usevardelta','var')
    usevardelta = 0;
end
if ~exist('usevartol','var')
    usevartol = 0;
end
    
if usevaralpha == 0 && usevarbeta == 0 && usevardelta == 0 && usevartol == 0
    close all;
    clearvars -except closed test inseti; % clear all;
    % closed: published code or closed
    % test  : which test used in tests
    % inseti: input of inseti-file
    
    usevaralpha = 0;
    usevarbeta = 0;
    usevardelta = 0;
    usevartol = 0;
    
    disp('Cleared all variables except of "inseti" in setInput.m');
end

%% Set parameters
%
%%
% *Description of "Set parameters"*
%
% Files with parameters are in the folder |inseti|.
% 
% The used parameter settings will be stored later in the folder |output|.
%
%%
% *Example of "Set parameters"*
%
%   inseti = 'example';
%   start;
%
%%
% *Input Arguments of "Set parameters"*
%
% * inseti  :   name of function including input parameters in folder inseti.
% * closed  :   closed is always 0 in public version 
%               (value 1 would mean that closed code is available)
%               (see also <init.html> for variable |closed|).
%
% *Output Arguments of "Set parameters"*
%
% * seti.inseti     :   name of function including input parameters in folder inseti.
% * seti.closed     :   is always 0 in public version 
%                       (value 1 would mean that closed code is available)
%                        (see also <init.html> for variable |closed|)
%
%%
% *Code*

if ~exist('inseti','var')
    inseti = '';
    % make sure to use '' instead of [] for empty variable in this case 
    % because eval(inseti) will not work with []
    fprintf('   Parameter inseti is empty.\n')
else
    fprintf('   Parameter inseti = %s.\n',inseti)
end

if ~exist('closed','var')
    closed = 0;
end
    
seti.inseti = inseti;
seti.closed = closed; % is set in init;
clear closed;

% make sure that a file exists and that the path to the file contains the string /inseti/.
% For Windows we look for the string \inseti\.
% Further check /insetiClosed/ and \insetiClosed\.
if isempty(inseti)
    fprintf('   Parameter inseti is not evaluated because it is empty.\n')
elseif exist(inseti,'file') == 2 && ...
        ( ~isempty(regexp(which(inseti),'/inseti/','once')) || ~isempty(regexp(which(inseti),'\inseti\','once')) ||...
          ~isempty(regexp(which(inseti),'/insetiClosed/','once')) || ~isempty(regexp(which(inseti),'\insetiClosed\','once')) )
    eval(inseti)
else
    fprintf('Error: %s is not a file in the folder /inseti/ or its subfolders.\n',inseti)
end

%% Various 2
%
%%
% *Description of "Various 2":
%
% Deal with various input parameters. Set the values for current run of start.
%
%%
% *Input Arguments of "Various 2"*
%
% Current values, which was set in <varalpha.html>, <varbeta.html>, <vardelta.html>
%
% * alphaVal    : alpha (in case of varalpha), see <varalpha.html>
% * betaVal     : beta (in case of varbeta), see <varbeta.html>
% * deltaVal    : delta (noise level) (in case of vardelta), see <vardelta.html>
% 
%%
% *See Also:*
%
% * <varalpha.html>
% * <varbeta.html>
% * <vardelta.html>
%%
% *Code*

if usevaralpha == 1
    seti.alpha = alphaVal;
end
if usevarbeta == 1
    seti.beta = betaVal;
end
if usevardelta == 1
    seti.delta = deltaVal;
end
if usevartol == 1
    seti.tol = tolVal;
end

%% Set and make directories for output
% make directories, see also <dirMake.html>.
%
%%
% *Input Arguments*
%
% All fields are set automatically, if they does not exist:
%
% * seti.dirOutput  :   name of folder in which directories for output 
%                       files and figures are created.
%                       (default: 'output')
% * seti.dirDatetime    :   date and time separated by the character "T",
%                           e.g. 20161006T105735.
% * seti.dirSuffix      :   suffix for dirname (default: '')
% * seti.dirSuffixAdd   :   additional suffix of dirname 
%                           (e.g. in case of varalpha, varbeta, vardelta, gscale; or expData is 'fresnel')
% * seti.dirname        :   dirname of files and figures for current
%                           computation, default:
%                           |seti.dirOutput + / + seti.datetime + _ + seti.inseti + seti.dirSuffix + seti.dirSuffixAdd|
%                           If seti.inseti is empty, 'noinseti' is used in
%                           dirname, e.g.
%                           |seti.dirname =
%                           'output/20161006T102740_noinseti'|.
% Not set automatically:
%
% * seti.expData : see <setData.html>
%
%%
% *Code*

seti = dirMake(seti,usevaralpha,usevarbeta,usevardelta,usevartol);

%% File suffix
% Set a suffix for files. In case of |varalpha|, |varbeta|, or |vardelta|
% it is done in these files.
%
% * seti.fileSuffix     :   Suffix for files (default: '').
%                           This is used in case of varalpha, varbeta, or vardelta.

seti = checkfield(seti,'fileSuffix','');

%% See Also
%
% * <start.html>
%
