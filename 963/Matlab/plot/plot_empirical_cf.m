cd('..');
boot;
cd(cgmm_config.directories.plot);

s = csvread(cgmm_config.time_series.file,1,0);
y = log(s);
r = diff(y);

% create an evaluation grid
grid_min = cgmm_config.cgmm.grid_min;
grid_res = cgmm_config.cgmm.grid_res;
grid_max = cgmm_config.cgmm.grid_max;
grid_margin = grid_min:(grid_max-grid_min)/(grid_res-1):grid_max;
omega = mgrid(grid_margin, 2); % evaluation grid

phi_emp = mean(exp(1i*r*omega)); % calculate empirical characteristic function

phi_emp = reshape(phi_emp, sqrt(length(phi_emp)), sqrt(length(phi_emp)));

subplot(1,2,1);
mesh(grid_margin, grid_margin, real(phi_emp));
title('Empirical - real part');
subplot(1,2,2);
mesh(grid_margin, grid_margin, imag(phi_emp));
title('Empirical - imaginary part');

print(cgmm_config.plots.cf('empirical'), cgmm_config.plots.device)
