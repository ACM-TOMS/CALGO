%% vardelta
% Call start with various delta
% (i.e. various input parameters for the noise level).
%
% Warning: This function clears variables, creates a directory and files.
%
%% Syntax
%
%   vardelta
%
%% Description 
%
% Call start with various input parameters for the noise level
% $\delta$ to tackle the inverse scattering problem with these parameters.
% This is useful to test easily with different noise levels $\delta$.
%
%% Example
%
% Set the values as a vector in this file, e.g.
%
%   delta = [10.^-4; 10.^-3];
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
% * <varbeta.html>
% * <vartol.html>
% * <start.html>
%
%% Code

disp(' ')
disp('----- Various delta (noise level)... -----')
disp(' ')

close all;
clearvars -except inseti;

%%
% *In-place setting of delta:*

% ----------------------------------------------------------------------
delta = [10.^-4; 10.^-3];
% ----------------------------------------------------------------------

%%
% *Process*

usevardelta = 1; %#ok (suppress MATLAB warning, because will be needed in start.m)
dirDatetime = datestr(now, 'yyyymmddTHHMMSS');
disDelta = zeros(length(delta),1);
errDelta = zeros(length(delta),1);
for i = 1:length(delta)
    seti.dirDatetime = dirDatetime;
    deltaVal = delta(i);
    % seed the random number generator
    seti.seed = floor(3800*pi*deltaVal);
    %
    seti.fileSuffix = sprintf('_iVar_%02d_delta_%g',i,deltaVal);
    %fprintf('%s:\n',seti.fileSuffix);
    start;
    close all;
    disDelta(i) = seti.dis(seti.iOutStop);
    errDelta(i) = seti.err(seti.iOutStop);
    clearvars -except inseti usevardelta dirDatetime delta i disDelta errDelta;
end
usevardelta = 0;

%%
% *Output*

disp(' ')
disp('Output')
disp('columns: iVar | delta | dis | err')
for i = 1:length(delta)
    fprintf('%02d | %g | %g | %g\n', i, delta(i), disDelta(i), errDelta(i))
end
