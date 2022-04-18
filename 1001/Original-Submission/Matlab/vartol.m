%% vartol
% Call start with various tolerance for GMRES.
% (i.e. various input parameters for the tolerance |seti.tol| tolerance for 
% GMRES in solveLippmannSchwinger.html, default: 1E-6).
%
% Warning: This function clears variables, creates a directory and files.
%
%% Syntax
%
%   vartol
%
%% Description 
%
% Call start with various input parameters for the tolerance |seti.tol| for 
% GMRES in <solveLippmannSchwinger.html> (default: 1E-6).
%
% This is useful to try different tolerances.
%
%% Example
%
% Set the values as a vector in this file, e.g.
%
%   tol = [1E-5; 1E-6; 1E-7];
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
% * <vardelta.html>
% * <start.html>
%
% Note that figures and results are stored in folder output.
%
%% Code

disp(' ')
disp('----- Various tol... -----')
disp(' ')

%%
close all;
clearvars -except inseti;

%%
% *In-place setting of tol:*

% ----------------------------------------------------------------------
tol = [1E-5; 1E-6; 1E-7];
% ----------------------------------------------------------------------

%%
% *Process*
%
% Set |usevartol| to 1 and store current |seti.tol| in |tolVal|. Call
% |start|.

usevartol = 1; %#ok (suppress MATLAB warning, because will be needed in start.m)
dirDatetime = datestr(now, 'yyyymmddTHHMMSS');
disTol = zeros(length(tol),1);
errTol = zeros(length(tol),1);
for i = 1:length(tol)
    seti.dirDatetime = dirDatetime;
    tolVal = tol(i);
    seti.fileSuffix = sprintf('_iVar_%02d_tol_%g',i,tolVal);
    %fprintf('%s:\n',seti.fileSuffix);
    start;
    close all;
    disTol(i) = seti.dis(seti.iOutStop);
    errTol(i) = seti.err(seti.iOutStop);
    clearvars -except inseti usevartol dirDatetime tol i disTol errTol;
end
usevartol = 0;

%%
% *Output*

disp(' ')
disp('Output')
disp('columns: iVar | tol | dis | err')
for i = 1:length(tol)
    fprintf('%02d | %g | %g | %g\n', i, tol(i), disTol(i), errTol(i))
end
