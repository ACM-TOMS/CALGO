cd('..');
boot;
cd(cgmm_config.directories.main);

% Init parameters
n = 2
p = 2
dt = 1/250;

% load PCSV parameters estimated from the corresponding time series
% to be used as 'true' parameters for simulation
load(cgmm_config.estimates.pcsv_partial);

% rename estimates to make clear that they play the role of true parameters here
mu = mu_cgmm
A = A_cgmm
lambda_0 = lambda_0_cgmm
kappa = kappa_cgmm
theta = theta_cgmm
sigma = sigma_cgmm
rho = rho_cgmm

% prepare simulation parameters
S_0 = [100 100];
y_0 = log(S_0);
time_steps = cgmm_config.monte_carlo.time_steps/dt;

% encode true parameters into flat parameter vector
[theta_flat_0, decode] = encode_pcsv_param_partial(mu, A, lambda_0, kappa ...
                                           , theta, sigma, rho);
% prepare characteristic function call with the encoded parameters
cf = @(omega, th, y_t, tau) cf_pcsv_partial_theta(decode, omega, th, y_t, tau);
% prepare options for optimization routine
options = optimset('Display', 'iter' ...
                  , 'Algorithm', 'interior-point');
% contraint for all parameters of +-10% around true parameters
lb = theta_flat_0 - abs(theta_flat_0)*0.1;
ub = theta_flat_0 + abs(theta_flat_0)*0.1;

simulation_runs = cgmm_config.monte_carlo.simulation_runs
if exist(cgmm_config.monte_carlo.pcsv_partial)
  load(cgmm_config.monte_carlo.pcsv_partial);
  % loaded file contains cell arrays 'first_step_estimates' and 'cgmm_estimates'
  first_run = size(first_step_estimates{1},1)
else
  first_step_estimates = {} % store a cell of first step estimates for each time step
  cgmm_estimates = {} % store a cell of cgmm estimates for each time step
  % each estimates{k} is a simulation_runs x number of params matrix
  for k = 1:length(time_steps)
    first_step_estimates{k} = theta_flat_0;
    cgmm_estimates{k} = theta_flat_0;
  end
  first_run = 1
end

for run = first_run:simulation_runs
  disp(strcat('Simulation run: ',num2str(run),'/',num2str(simulation_runs)));
  t = 0:dt:(time_steps(end)*dt);
  % simulate time series
  [y, lambda] = sim_pcsv_partial(y_0, mu, A, lambda_0, kappa, theta, sigma, rho, t);
  for k = 1:length(time_steps)
    y_k = y(1:time_steps(k),:);
    % perform parameter estimation
    tic;
    [theta_flat_cgmm, theta_flat_first] = cgmm(y_k, dt, cf, theta_flat_0 ...
                                            , cgmm_config.cgmm.grid_min+1 ...
                                            , cgmm_config.cgmm.grid_max+1 ...
                                            , cgmm_config.cgmm.grid_res ...
                                            , lb, ub, options);
    toc;
    first_step_estimates{k} = [first_step_estimates{k}; theta_flat_first];
    cgmm_estimates{k} = [cgmm_estimates{k}; theta_flat_cgmm];
  end
  % save monte carlo estimates in each iteration
  save(cgmm_config.monte_carlo.pcsv_partial, 'first_step_estimates' ...
       , 'cgmm_estimates', 'time_steps');
end
