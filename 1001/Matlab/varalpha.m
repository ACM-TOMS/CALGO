%% varalpha
% Call start with various alpha
% (i.e. various input parameters for the regularization parameter alpha).
%
% Warning: This function clears variables, creates a directory and files.
%
%% Syntax
%
%   varalpha
%
%% Description 
%
% Call start with various input parameters for the regularization parameter
% $\alpha$ (sparsity) to tackle the inverse scattering problem with these parameters.
% This is useful to find a suitable $\alpha$.
%
%% Example
%
% Set the values as a vector in this file, e.g.
%
%   alpha = [100; 500; 1000];
%
% In terminal:
%
%   inseti = 'example'; % setting of seti.alpha inside is overwrittten by varalpha
%   varalpha
%
% See also <start.html> in section "Advanced starts".
%
%% Input Arguments
%
% * inseti    : name of file with input parameters in folder |inseti|.
%
% Set |inseti| before. Note that all other parameters will be cleared.
%
%% In-place Input
% Define in this file.
%
%   alpha     : vector with values of several regularization parameters
%               (the current parameter for computation is stored in
%               |alphaVal|)
%
%% Output Arguments
%
% * usevaralpha                 : 0 
%                                (Inside this file it is 1, i.e. we
%                                want to use various alpha inputs,
%                                afterwars it is 0.)
% * dirDatetime and seti.dirDatetime  : date and time in format yyyymmddTHHMMSS, e.g. 20160815T105311
% * errAlpha                    : vector with relative errors of reconstruction depending on used alpha
%
% Note that figures and results are stored in folder output.
%
%% More About
%
% <varalpha.html>, <varbeta.html>, <vardelta.html> and <vartol.html>
% have the same structure:
%
% * alphaVal    : alpha (in case of varalpha), see <varalpha.html>
% * betaVal     : beta (in case of varbeta), see <varbeta.html>
% * deltaVal    : delta (noise level) (in case of vardelta), see <vardelta.html>
% * tolVal      : tol (tolerance of GMRES) (in case of vartol), see <vartol.html>.
% 
% Vectors containing relative reconstruction errors...
%
% * errAlpha    : ... for each alpha (in case of varalpha).
% * errBeta     : ... for each beta (in case of varbeta).
% * errDelta    : ... for each delta (noise level) (in case of vardelta).
% * errTol      : ... for each tolerance for GMRES (in case of vartol).
%
%% See Also
%
% * <setInput.html>
% * <varbeta.html>
% * <vardelta.html>
% * <vartol.html>
%
%% Code

disp(' ')
disp('----- Various alpha... -----')
disp(' ')

%%
close all;
clearvars -except inseti;

%%
% *In-place setting of alpha:*

% ----------------------------------------------------------------------
alpha = [100; 500; 1000];
% ----------------------------------------------------------------------

%%
% *Process*
%
% Set |usevaralpha| to 1 and store current $\alpha$ in |alphaVal|. Call
% |start|.

usevaralpha = 1; %#ok (suppress MATLAB warning, because will be needed in start.m)
dirDatetime = datestr(now, 'yyyymmddTHHMMSS');
disAlpha = zeros(length(alpha),1);
errAlpha = zeros(length(alpha),1);
for i = 1:length(alpha)
    seti.dirDatetime = dirDatetime;
    alphaVal = alpha(i);
    seti.fileSuffix = sprintf('_iVar_%02d_alpha_%g',i,alphaVal);
    %fprintf('%s:\n',seti.fileSuffix);
    start;
    close all;
    disAlpha(i) = seti.dis(seti.iOutStop);
    errAlpha(i) = seti.err(seti.iOutStop);
    clearvars -except inseti usevaralpha dirDatetime alpha i disAlpha errAlpha;
end
usevaralpha = 0;

%%
% *Output*

disp(' ')
disp('Output')
disp('columns: iVar | alpha | dis | err')
for i = 1:length(alpha)
    fprintf('%02d | %g | %g | %g\n', i, alpha(i), disAlpha(i), errAlpha(i))
end
