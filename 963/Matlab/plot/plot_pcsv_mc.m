cd('..');
boot;
cd(cgmm_config.directories.plot);

load(cgmm_config.monte_carlo.pcsv);

all_re_cgmm = zeros(size(cgmm_estimates{1},1)-1,length(time_steps));
all_re_fs = zeros(size(cgmm_estimates{1},1)-1,length(time_steps));

mean_mre_cgmm = zeros(1, length(time_steps));
mean_mre_fs = zeros(1, length(time_steps));
std_mre_cgmm = zeros(1, length(time_steps));
std_mre_fs = zeros(1, length(time_steps));

for k=1:length(time_steps)
  theta_fs = first_step_estimates{k};
  theta_true = theta_fs(1,:);
  theta = theta_fs(2:end,:);
  
  theta_true = repmat(theta_true, size(theta,1) ,1);
  re = abs(theta-theta_true) ./ abs(theta_true);
  mre = mean(re);
  stdre = std(re);
  % all estimates
  all_re_fs(:,k) = mean(re,2);
  mean_mre_fs(k) = mean(mre);
  std_mre_fs(k) = std(mre);
  % kappa
  kappa_mre_fs(k) = mean(mre(1:2));
  % theta
  theta_mre_fs(k) = mean(mre(3:4));
  % sigma
  sigma_mre_fs(k) = mean(mre(5:6));
  % rho
  rho_mre_fs(k) = mean(mre(7:8));
  
  theta_cgmm = cgmm_estimates{k};
  theta = theta_cgmm(2:end,:);
  re = abs(theta-theta_true) ./ abs(theta_true);
  mre = mean(re);
  stdre = std(re);
  % all estimates
  all_re_cgmm(:,k) = mean(re,2);
  mean_mre_cgmm(k) = mean(mre);
  std_mre_cgmm(k) = std(mre);
  % kappa
  kappa_mre_cgmm(k) = mean(mre(1:2));
  % theta
  theta_mre_cgmm(k) = mean(mre(3:4));
  % sigma
  sigma_mre_cgmm(k) = mean(mre(5:6));
  % rho
  rho_mre_cgmm(k) = mean(mre(7:8));
end

all_re_cgmm(:,end) = NA;
mean_mre_cgmm(end) = NA;
std_mre_cgmm(end) = NA;
% kappa
kappa_mre_cgmm(end) = NA;
% theta
theta_mre_cgmm(end) = NA;
% sigma
sigma_mre_cgmm(end) = NA;
% rho
rho_mre_cgmm(end) = NA;

subplot(2,1,1);

plot(time_steps, all_re_fs, 'k');
hold on;
plot(time_steps, mean_mre_fs, 'r', 'linewidth', 5);
hold off;
legend('100 independent simulations', 'Mean of all simulations', 'location', 'southwest');
ylabel('mean relative error');
xlabel('number of data points in time series');
title('PCSV CMM mean relative estimation errors');

subplot(2,1,2);
plot(time_steps, std_mre_fs);
ylabel('standard deviation of mean relative error');
xlabel('number of data points in time series');

print(cgmm_config.plots.mre.all('pcsv'), cgmm_config.plots.device)

subplot(2,2,1);
plot(time_steps, kappa_mre_fs);
title('kappa');
subplot(2,2,2);
plot(time_steps, theta_mre_fs);
title('theta');
subplot(2,2,3);
plot(time_steps, sigma_mre_fs);
title('sigma');
subplot(2,2,4);
plot(time_steps, rho_mre_fs);
title('rho');

print(cgmm_config.plots.mre.individual('pcsv'), cgmm_config.plots.device)
