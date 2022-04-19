cd('..');
boot;
cd(cgmm_config.directories.plot);

load(cgmm_config.monte_carlo.wasc);

all_re = zeros(size(estimates{1},1)-1,length(time_steps));

mean_mre = zeros(1, length(time_steps));
std_mre = zeros(1, length(time_steps));

for k=1:length(time_steps)
  theta = estimates{k};
  theta_true = theta(1,:);
  theta = theta(2:end,:);
  theta_true = repmat(theta_true, size(theta,1) ,1);
  re = abs(theta-theta_true) ./ abs(theta_true);
  mre = mean(re);
  stdre = std(re);
  % all estimates
  all_re(:,k) = mean(re,2);
  mean_mre(k) = mean(mre);
  std_mre(k) = std(mre);
  % M
  M_mre(k) = mean(mre(1:4));
  % Q
  Q_mre(k) = mean(mre(5:8));
  % rho
  rho_mre(k) = mean(mre(9:10));
  % beta
  beta_mre(k) = mre(11);
end

subplot(2,1,1);

plot(time_steps, all_re, 'k');
hold on;
plot(time_steps, mean_mre, 'r', 'linewidth', 5);
hold off;
legend('100 independent simulations', 'Mean of all simulations', 'location', 'southwest');
ylabel('mean relative error');
xlabel('number of data points in time series');
title('WASC first step mean relative estimation errors');

subplot(2,1,2);
plot(time_steps, std_mre);
ylabel('standard deviation of mean relative error');
xlabel('number of data points in time series');

print(cgmm_config.plots.mre.all('wasc'), cgmm_config.plots.device)

subplot(2,2,1);
plot(time_steps, M_mre);
title('M');
subplot(2,2,2);
plot(time_steps, Q_mre);
title('Q');
subplot(2,2,3);
plot(time_steps, rho_mre);
title('rho');
subplot(2,2,4);
plot(time_steps, beta_mre);
title('beta');

print(cgmm_config.plots.mre.individual('wasc'), cgmm_config.plots.device)
