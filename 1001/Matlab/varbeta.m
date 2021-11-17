%% varbeta
% Call start with various beta
% (i.e. various input parameters for the regularization parameter beta).
%
% Warning: This function clears variables, creates a directory and files.
%
%% Syntax
%
%   varbeta
%
%% Description 
%
% Call start with various input parameters for the regularization parameter
% $\beta$ (total variation) to tackle the inverse scattering problem with these parameters.
% This is useful to find a suitable $\beta$.
%
%% Example
%
% Set the values as a vector in this file, e.g.
%
%   beta = [6E-5; 8E-5; 2E-4; 3E-4; 4E-4];
%
% See <varalpha.html> or <start.html> for usage.
%
%% Input and Output Arguments
%
% Input and Ouptut is analogous to <varalpha.html>.
%
%% See Also
%
% * <setInput.html>
% * <varalpha.html>
% * <vardelta.html>
% * <vartol.html>
% * <start.html>
%
%% Code

disp(' ')
disp('----- Various beta... -----')
disp(' ')

close all;
clearvars -except inseti;

%%
% *In-place setting of beta:*

% ----------------------------------------------------------------------
beta = [6E-5; 8E-5; 2E-4; 3E-4; 4E-4];
% ----------------------------------------------------------------------

%%
% *Process*

usevarbeta = 1; %#ok (suppress MATLAB warning, because will be needed in start.m)
dirDatetime = datestr(now, 'yyyymmddTHHMMSS');
disBeta = zeros(length(beta),1);
errBeta = zeros(length(beta),1);
for i = 1:length(beta)
    seti.dirDatetime = dirDatetime;
    betaVal = beta(i);
    seti.fileSuffix = sprintf('_iVar_%02d_beta_%g',i,betaVal);
    %fprintf('%s:\n',seti.fileSuffix);
    start;
    close all;
    disBeta(i) = seti.dis(seti.iOutStop);
    errBeta(i) = seti.err(seti.iOutStop);
    clearvars -except inseti usevarbeta dirDatetime beta i disBeta errBeta;
end
usevarbeta = 0;

%%
% *Output*

disp(' ')
disp('Output')
disp('columns: iVar | beta | dis | err')
for i = 1:length(beta)
    fprintf('%02d | %g | %g | %g\n', i, beta(i), disBeta(i), errBeta(i))
end
