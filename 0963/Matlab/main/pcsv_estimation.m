cd('..');
boot;
cd(cgmm_config.directories.main);

% Init parameters
n = 2
p = 2
dt = 1/250;

% load two-dimensional time series for estimation
s = csvread(cgmm_config.time_series.file,1,0);
y = log(s);
r_real = diff(y);
y_0 = y(1,:);

% perform heuristic estimation of parameters
[mu_heur, A_heur, lambda_0_heur, kappa_heur, theta_heur ...
  , sigma_heur, rho_heur] = heuristic_pcsv_param(y, dt, p)

% encode parameters into flat parameter vector
% call it theta_flat to avoid naming conflicts with the mean reversion
% level parameter theta
[theta_flat_0, decode] = encode_pcsv_param(mu_heur, A_heur, lambda_0_heur ...
  , kappa_heur, theta_heur, sigma_heur, rho_heur);

% prepare characteristic function call with the encoded parameters
cf = @(omega, th, y_t, tau) cf_pcsv_v_theta(decode, omega, th, y_t, tau);

% prepare options for optimization routine
options = optimset('Display', 'iter' ...
                  , 'Algorithm', 'interior-point');
% contraint for all parameters of +-25% around heuristic estimate
lb = theta_flat_0 - abs(theta_flat_0)*0.25;
ub = theta_flat_0 + abs(theta_flat_0)*0.25;
% except for correlation (where absolute constraints are more feasible)
lb(end-1:end) = max(theta_flat_0(end-1:end)-0.6, -1);
ub(end-1:end) = min(theta_flat_0(end-1:end)+0.6, 1);

% test for the Feller condition
feller_condition = @(kappa, theta, sigma) 2*kappa.*theta > sigma.^2;
if ~all(feller_condition(kappa_heur, theta_heur, sigma_heur))
  error('Heuristic estimates do not satisfy the Feller condition!');
end

% set parameter constraints for kappa, theta and sigma
% so that the Feller condition will still be satisfied after optimization:
% -> set lower bounds for kappa and theta
% -> set upper bounds for sigma
% distribute constraints "relatively equal" among the parameters, i.e.
% (1-alpha)^2 * 2 * kappa * theta > (1+alpha)^2 * sigma^2
% If the Feller condition holds, alpha is positive.
alpha = ( sqrt(2*kappa_heur.*theta_heur) - sigma_heur ) ./ ...
        ( sqrt(2*kappa_heur.*theta_heur) + sigma_heur );

idx_kappa = 1:2;
idx_theta = 3:4;
idx_sigma = 5:6;

lb(idx_kappa) = max(lb(idx_kappa), (theta_flat_0(idx_kappa)-abs(theta_flat_0(idx_kappa)).*alpha));
lb(idx_theta) = max(lb(idx_theta), (theta_flat_0(idx_theta)-abs(theta_flat_0(idx_theta)).*alpha));
ub(idx_sigma) = min(ub(idx_sigma), (theta_flat_0(idx_sigma)+abs(theta_flat_0(idx_sigma)).*alpha));

% perform parameter estimation
tic;
[theta_flat_cgmm, theta_flat_first] = cgmm(y, dt, cf, theta_flat_0 ...
                                            , cgmm_config.cgmm.grid_min+1 ...
                                            , cgmm_config.cgmm.grid_max+1 ...
                                            , cgmm_config.cgmm.grid_res ...
                                            , lb, ub, options);
toc;

% decode parameters
[mu_first, A_first, lambda_0_first, kappa_first, theta_first ...
  , sigma_first, rho_first] = decode_pcsv_param(theta_flat_first, decode);

[mu_cgmm, A_cgmm, lambda_0_cgmm, kappa_cgmm, theta_cgmm ...
  , sigma_cgmm, rho_cgmm] = decode_pcsv_param(theta_flat_cgmm, decode);

% save heuristic estimates, first step estimates and cgmm estimates
save( ...
  cgmm_config.estimates.pcsv ...
  , 'theta_flat_0', 'mu_heur', 'A_heur', 'lambda_0_heur' ...
  , 'kappa_heur', 'theta_heur', 'sigma_heur', 'rho_heur' ...
  , 'theta_flat_first', 'mu_first', 'A_first', 'lambda_0_first' ...
  , 'kappa_first', 'theta_first', 'sigma_first', 'rho_first' ...
  , 'theta_flat_cgmm', 'mu_cgmm', 'A_cgmm', 'lambda_0_cgmm' ...
  , 'kappa_cgmm', 'theta_cgmm', 'sigma_cgmm', 'rho_cgmm' ...
)
