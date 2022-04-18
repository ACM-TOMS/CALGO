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

S_0 = repmat(100,1,n);
y_0 = log(S_0);

% create an evaluation grid
grid_min = cgmm_config.cgmm.grid_min;
grid_res = cgmm_config.cgmm.grid_res;
grid_max = cgmm_config.cgmm.grid_max;
margin = grid_min:(grid_max-grid_min)/(grid_res-1):grid_max;
omega = mgrid(margin, 2); % evaluation grid

% simulate
dt = 1/250
N = 500;
phi_emp = cf_pcsv_sim(mu, A, lambda_0, kappa, theta, sigma, rho, omega, repmat(y_0, 1/dt, 1), dt, N);
phi_emp = mean(phi_emp);

phi = cf_pcsv_alternative(mu, A, lambda_0, kappa, theta, sigma, rho, omega, repmat(y_0, 1/dt, 1), dt);
phi = mean(phi);

phi = reshape(phi, sqrt(length(phi)), sqrt(length(phi)));
phi_emp = reshape(phi_emp, sqrt(length(phi_emp)), sqrt(length(phi_emp)));

subplot(2,2,1);
mesh(margin, margin, real(phi));
h = title('Theoretical - real part');

subplot(2,2,3);
mesh(margin, margin, real(phi_emp));
title('Empirical - real part');
subplot(2,2,2);
mesh(margin, margin, imag(phi));
title('Theoretical - imaginary part');
subplot(2,2,4);
mesh(margin, margin, imag(phi_emp));
title('Empirical - imaginary part');

print('../output/png/pcsv_cf.png', cgmm_config.plots.device)
