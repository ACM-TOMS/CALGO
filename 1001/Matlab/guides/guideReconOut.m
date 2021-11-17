% This example does the same as guideReconDetails.m, 
% but plots and saves figures, further saves files.

inseti = 'exampleMod';
init;                       % Initialization (addpath...)
setInput;                   % Close figures, clear variables, create folder for output.
%ticGuide = tic;
seti = setData(seti,4,2);   % Set data: experimental set-up, contrast, simulate data.
seti = setRecon(seti,4,2);  % Settings for variational reconstruction.
seti = recon(seti,4,2);     % Variational reconstruction (process).
%tocGuide = toc(ticGuide)
