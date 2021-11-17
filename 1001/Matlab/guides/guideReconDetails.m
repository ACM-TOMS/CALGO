% This example does the same as guideRecon.m, but displays more messages.

init;                     % Initialization (addpath...)
seti = struct;
seti = setData(seti,4);   % Set data: experimental set-up, contrast, simulate data
seti = setRecon(seti,4);  % Settings for variational reconstruction
seti = recon(seti,4);     % Variational reconstruction (process)
