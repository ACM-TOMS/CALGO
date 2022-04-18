cd('..');
boot;
cd(cgmm_config.directories.plot);

load(cgmm_config.estimates.pcsv);

n = 2
p = 2
dt = 1/250

% rename estimates to make clear that they play the role of true parameters here
mu = mu_cgmm
A = A_cgmm
lambda_0 = lambda_0_cgmm
kappa = kappa_cgmm
theta = theta_cgmm
sigma = sigma_cgmm
rho = rho_cgmm

s = csvread(cgmm_config.time_series.file,1,1);
y = log(s);
r = diff(y);

% create an evaluation grid
grid_min = cgmm_config.cgmm.grid_min;
grid_res = cgmm_config.cgmm.grid_res;
grid_max = cgmm_config.cgmm.grid_max;
grid_margin = grid_min:(grid_max-grid_min)/(grid_res-1):grid_max;
omega = mgrid(grid_margin, 2); % evaluation grid

phi_emp = mean(exp(1i*r*omega)); % calculate empirical characteristic function
phi = cf_pcsv_v(mu, A, lambda_0, kappa, theta, sigma, rho, omega, zeros(30,1), dt);
phi = phi(end,:);

phi = reshape(phi, sqrt(length(phi)), sqrt(length(phi)));
phi_emp = reshape(phi_emp, sqrt(length(phi_emp)), sqrt(length(phi_emp)));

subplot(2,2,1);
mesh(grid_margin, grid_margin, real(phi));
h = title('Theoretical - real part');
P = get(h,'Position');
set(h,'Position',[P(1)+50 P(2)+100 P(3)+40])
set(h,'visible', 'off');

subplot(2,2,3);
mesh(grid_margin, grid_margin, real(phi_emp));
title('Empirical - real part');
subplot(2,2,2);
mesh(grid_margin, grid_margin, imag(phi));
title('Theoretical - imaginary part');
subplot(2,2,4);
mesh(grid_margin, grid_margin, imag(phi_emp));
title('Empirical - imaginary part');

print(cgmm_config.plots.cf('pcsv'), cgmm_config.plots.device)
