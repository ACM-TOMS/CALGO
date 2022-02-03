cd('..');
boot;
cd(cgmm_config.directories.plot);

load(cgmm_config.estimates.wasc);

% rename estimates to make clear that they play the role of true parameters here
mu = mu_cgmm
Sigma_0 = Sigma_0_cgmm
M = M_cgmm
Q = Q_cgmm
rho = rho_cgmm
beta = round(beta_cgmm) % integer beta required for simulation

S_0 = [100; 100]
y_0 = log(S_0)

dt = 1/250
t = 0:dt:5;

[y, vol] = sim_wasc_2d(y_0, mu, Sigma_0, M, Q, rho, beta, t);

subplot(3,1,1);
plot(exp(y));
title('Stock Price');
subplot(3,1,2);
plot(vol(:,1:2));
title('Volatility');
subplot(3,1,3);
plot(vol(:,3));
title('Correlation');

print(cgmm_config.plots.simulation('wasc'), cgmm_config.plots.device)
