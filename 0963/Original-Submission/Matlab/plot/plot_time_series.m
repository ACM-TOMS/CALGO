cd('..');
boot;
cd(cgmm_config.directories.plot);

% load two-dimensional time series for estimation
s = csvread(cgmm_config.time_series.file,1,0);

plot(s)

print([cgmm_config.directories.png '/time_series' ...
    cgmm_config.time_series.name '.png'], cgmm_config.plots.device)
