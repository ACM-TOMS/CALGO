cd('..');
boot;
cd(cgmm_config.directories.main);

% Init parameters
n = 2
dt = 1/250;

% load two-dimensional time series for estimation
s = csvread(cgmm_config.time_series.file,1,0);
y = log(s);
r_real = diff(y);
y_0 = y(1,:);

% perform heuristic estimation of parameters
[mu_heur, Sigma_0_heur, M_heur, Q_heur, rho_heur, beta_heur] ...
      = heuristic_wasc_param_2d(y, dt)

% encode parameters into flat parameter vector
% call it theta_flat to avoid naming conflicts with the mean reversion
% level parameter theta
[theta_flat_0, decode] = encode_wasc_param(mu_heur, Sigma_0_heur, M_heur ...
                                  , Q_heur, rho_heur, beta_heur);

% prepare characteristic function call with the encoded parameters
cf = @(omega, th, y_t, tau) cf_wasc_v_theta(decode, omega, th, y_t, tau);

% prepare options for optimization routine
options = optimset('Display', 'iter' ...
                  , 'Algorithm', 'interior-point');
% contraint for all parameters of +-25% around heuristic estimate
lb = theta_flat_0 - abs(theta_flat_0)*0.25;
ub = theta_flat_0 + abs(theta_flat_0)*0.25;
lb(6) = -max(abs(theta_flat_0(5:8)))*0.25; % lower left Q; heurist estimate = 0
ub(6) = max(abs(theta_flat_0(5:8)))*0.25; % lower left Q; heurist estimate = 0
% except for correlation (where absolute constraints are more feasible)
lb(end-2:end-1) = max(theta_flat_0(end-2:end-1)-0.6, -1);
ub(end-2:end-1) = min(theta_flat_0(end-2:end-1)+0.6, 1);
% and beta
lb(end) = 7;
ub(end) = 13;

% perform parameter estimation
tic;
[theta_flat_cgmm, theta_flat_first] = cgmm(y, dt, cf, theta_flat_0 ...
                                            , cgmm_config.cgmm.grid_min ...
                                            , cgmm_config.cgmm.grid_max ...
                                            , cgmm_config.cgmm.grid_res ...
                                            , lb, ub, options);
toc;

% decode parameters
[mu_first, Sigma_0_first, M_first, Q_first, rho_first, beta_first] ...
                = decode_wasc_param(theta_flat_first, decode);

[mu_cgmm, Sigma_0_cgmm, M_cgmm, Q_cgmm, rho_cgmm, beta_cgmm] ...
                = decode_wasc_param(theta_flat_cgmm, decode);

% save heuristic estimates, first step estimates and cgmm estimates
save( ...
  cgmm_config.estimates.wasc ...
  , 'theta_flat_0', 'mu_heur', 'Sigma_0_heur', 'M_heur' ...
  , 'Q_heur', 'rho_heur', 'beta_heur' ...
  , 'theta_flat_first', 'mu_first', 'Sigma_0_first', 'M_first' ...
  , 'Q_first', 'rho_first', 'beta_first' ...
  , 'theta_flat_cgmm', 'mu_cgmm', 'Sigma_0_cgmm', 'M_cgmm' ...
  , 'Q_cgmm', 'rho_cgmm', 'beta_cgmm' ...
)
