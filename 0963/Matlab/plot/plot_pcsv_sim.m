cd('..');
boot;
cd(cgmm_config.directories.plot);

load(cgmm_config.estimates.pcsv);

n = 2
p = 2

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
kappa = 2*kappa
disp('Feller condition');
2*kappa.*theta >= sigma.^2

t = 0:1/250:1;
[y, lambda] = sim_pcsv(y_0, mu, A, lambda_0, kappa, theta, sigma, rho, t);

v = zeros(size(lambda,1),2);
c = zeros(size(lambda,1),1);
for k=1:length(c)
  M = A(:,1:p)*diag(lambda(k,:))*A(:,1:p)';
  v(k,:) = [sqrt(M(1,1)) sqrt(M(2,2))];
  c(k) = M(2,1) / (sqrt(M(1,1)*M(2,2)));
end

% plot simulated paths
subplot(3,1,1);
plot(exp(y));
title('Stock Price');
subplot(3,1,2);
plot(v);
title('Volatility');
subplot(3,1,3);
plot(c);
title('Correlation');

%print(cgmm_config.plots.feature_test('pcsv'), cgmm_config.plots.device)
