cd('..');
boot;
cd(cgmm_config.directories.main);

% Init parameters
n = 2
dt = 1/250;

% load WASC parameters estimated from the corresponding time series
% to be used as 'true' parameters for simulation
load(cgmm_config.estimates.wasc);

% rename estimates to make clear that they play the role of true parameters here
mu = mu_cgmm
Sigma_0 = Sigma_0_cgmm
M = M_cgmm
Q = Q_cgmm
rho = rho_cgmm
beta = round(beta_cgmm) % integer beta required for simulation

% prepare simulation parameters
S_0 = [100 100];
y_0 = log(S_0);
time_steps = cgmm_config.monte_carlo.time_steps/dt;

% encode true parameters into flat parameter vector
[theta_flat_0, decode] = encode_wasc_param(mu, Sigma_0, M, Q, rho, beta);
% prepare characteristic function call with the encoded parameters
cf = @(omega, th, y_t, tau) cf_wasc_v_theta(decode, omega, th, y_t, tau);
% prepare options for optimization routine
options = optimset('Display', 'iter' ...
                  , 'Algorithm', 'interior-point');
% contraint for all parameters of +-10% around heuristic estimate
lb = theta_flat_0 - abs(theta_flat_0)*0.1;
ub = theta_flat_0 + abs(theta_flat_0)*0.1;
%lb(6) = -max(abs(theta_flat_0(5:8)))*0.1; % lower left Q; heurist estimate = 0
%ub(6) = max(abs(theta_flat_0(5:8)))*0.1; % lower left Q; heurist estimate = 0
% except for correlation (where absolute constraints are more feasible)
%lb(end-2:end-1) = max(theta_flat_0(end-2:end-1)-0.6, -1);
%ub(end-2:end-1) = min(theta_flat_0(end-2:end-1)+0.6, 1);

simulation_runs = cgmm_config.monte_carlo.simulation_runs
if exist(cgmm_config.monte_carlo.wasc)
  load(cgmm_config.monte_carlo.wasc);
  % loaded file contains cell array 'estimates'
  first_run = size(estimates{1},1) + 1
else
  estimates = {}; % store a cell of estimates for each time step
  % each estimates{k} is a simulation_runs x number of params matrix
  for k = 1:length(time_steps)
    estimates{k} = theta_flat_0;
  end
  first_run = 1
end

for run = first_run:simulation_runs
  disp(strcat('Simulation run: ',num2str(run),'/',num2str(simulation_runs)));
  for k = 1:length(time_steps)
    t = 0:dt:(time_steps(k)*dt);
    % simulate time series
    y = sim_wasc_2d(y_0, mu, Sigma_0, M, Q, rho, beta, t);
    % perform parameter estimation
    tic;
    [theta_flat_cgmm, theta_flat_first] = cgmm(y, dt, cf, theta_flat_0 ...
                                            , cgmm_config.cgmm.grid_min ...
                                            , cgmm_config.cgmm.grid_max ...
                                            , cgmm_config.cgmm.grid_res ...
                                            , lb, ub, options);
    %theta_flat_cgmm = theta_flat_first = theta_flat_0;
    toc;
    estimates{k} = [estimates{k}; theta_flat_first];
  end
  % save monte carlo estimates in each iteration
  save(cgmm_config.monte_carlo.wasc, 'estimates', 'time_steps');
end

