% this file replaces file with same name at \core\estimation

function [theta_cgmm, theta_first] = cgmm(y, tau, cf, theta_0 ...
    , grid_min, grid_max, grid_res, lb, ub, options)
%CGMM Performs continuum of moments estimation for the model specified by cf.
%
%  [theta_cgmm, theta_first] = cgmm(y, tau, cf, theta_0, grid_min, grid_max
%    , grid_res, lb, ub, options) performs conditional continuum of moments
%    estimation for the model specified by its characteristic function cf,
%    based on the time series of log prices y and a start parameter theta_0
%    for the optimization.
%
%    INPUT       y: log prices of the assets
%              tau: real valued time difference
%               cf: characteristic function. cf is supposed to expect omega, 
%                   theta, y_t, tau with
%                   - omega a size(y,2) x grid_res matrix
%                   - theta the flat parameter vector (interpretation by cf)
%                   - y_t a (size(y,1)-1) x size(y,2) matrix
%                   - tau a real valued time difference
%                   cf is supposed to calculate to characteristic function
%                   of the continuous returns conditional on the log prices.
%          theta_0: initial parameter vector
%         grid_min: left end of marginal evaluation grid
%         grid_max: right end of marginal evaluation grid
%         grid_res: resolution of marginal evaluation grid
%               lb: lower bounds for the parameters
%               ub: upper bounds for the parameters
%          options: structure (or vector) of options for the optimization
%
%    OUTPUT  theta_cgmm: CGMM parameter estimate
%           theta_first: first step parameter estimate
%
% See also CF_PCSV_V_THETA, CF_WASC_THETA.
%
% created by Benedikt Rudolph
% DATE: 20-Aug-2012
  
  p = size(y,2); % dimension of the model
  T = size(y,1); % number of observations
  
  % choose a fixed grid according to the model's dimension for integration
  grid_margin = grid_min:(grid_max-grid_min)/(grid_res-1):grid_max;
  omega = mgrid(grid_margin, p);
  
  % prepare fixed parameters
  r = diff(y);
  s = (grid_max-grid_min)/6;
  % mv std normal
  w_kernel = @(x) ones(size(x,2),1);% exp(-0.5*diag(x'*x)/(s^2))/((2*3.1415*s)^(p*0.5));
  pi = w_kernel(omega)'; % weighting for L^2(pi) norm
  pi = pi/sum(pi);
  sqrt_pi = sqrt(pi);
  
  % precalculate empirical characteristic function (independent from theta)
  ecf = exp(1i*r*omega);
  
  % obtain the theoretical characteristic function as a function of theta
  phi = @(theta) cf(omega, theta, y, tau);
  
  % define the moment function depending on the parameters (theta)
  H = @(theta) ecf - phi(theta);
  
  % define the error function for the first step estimator
  % i.e. the I_{T-1} matrix as weighting matrix
  first_step_error = @(theta) mean(mean(H(theta)).*conj(mean(H(theta)))); 
  
  % perform first step optimization
  disp('Calculating first step estimator...');
  first_step_opt = prepare_nl_opt_call(first_step_error, theta_0 ...
                                      , lb, ub, options);
  theta_first = first_step_opt();
  theta_first = reshape(theta_first, 1, length(theta_first));
  disp('done');
  
  % CGMM estimation for correlated moment functions
  disp('Calculating optimal weighting matrix...');
  H_first = H(theta_first);
  H_first_weighted = H_first.*repmat(sqrt_pi,T-1,1);
  S = T/6;
  ph_kernel = @(x) ones(size(x,2),1);

  ph = ph_kernel((0:T)/S);
  ph = ph/sum(ph);
  UH_first = ph(1)*conj(H_first);
  ph = repmat(ph, 1, size(H_first,2));
  for t=1:(T-1)
    UH_first(t,:) = UH_first(t,:) ...
                    + sum(ph(2:t,:).*conj(H_first(1:(t-1),:)), 1) ...
                    + sum(ph(2:(T-t),:).*conj(H_first((t+1):(T-1),:)), 1);
  end
  UH_first_weighted = UH_first.*repmat(sqrt_pi,T-1,1);
  C = (UH_first_weighted*H_first_weighted')/(T-length(theta_first));
  alpha = 0.02; % regularization parameter
  opt_weighting = inv(alpha*eye(T-1)+C^2);
  disp('done');
  w = @(mean_H_theta) H_first_weighted * (mean_H_theta.*sqrt_pi)';
  v = @(mean_H_theta) UH_first_weighted * (mean_H_theta.*sqrt_pi).';
  % define the second step error function for the efficient estimator
  % wrap the call to mean(H(theta)) to prevent calling it twice (expensive!)
  wrap_second_step_error = @(mean_H_theta) real(w(mean_H_theta)' ...
                      *opt_weighting*v(mean_H_theta));
  second_step_error = @(theta) wrap_second_step_error(mean(H(theta)));
  
  % prepare CGMM optimization
  second_step_opt = prepare_nl_opt_call(second_step_error, theta_first ...
                                      , lb, ub, options);
  
  disp('Calculating second step cgmm estimator...');
  % call to optimization function
  theta_cgmm = second_step_opt();
  
  theta_cgmm = reshape(theta_cgmm, 1, length(theta_cgmm));
  disp('done');
end
