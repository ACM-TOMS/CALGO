function [y, lambda] = sim_pcsv(y_0, mu, A, lambda_0, kappa ...
                               , theta, sigma, rho, t)
%SIM_PCSV Simulates a PCSV model for the given time points in t.
%
%  [y, lambda] = sim_pcsv(y_0, mu, A, lambda_0, kappa, theta, sigma, rho, t)
%    Performs an Euler simulation of the PCSV model specified by the given
%    parameters for the time points in t.
%    Assuming that n is the dimension of the problem, p the number of
%    eigenvalues/eigenvectors, the parameters expected are
%
%    INPUT      y_0: n-dimensional vector with initial log prices
%                mu: n-dimensional mean vector
%                 A: n x p matrix of eigenvectors
%          lambda_0: p-dimensional vector of initial eigenvalues
%             kappa: p-dimensional vector representing "mean reversion speed"
%             theta: p-dimensional vector of "mean reversion means"
%             sigma: p-dimensional vector of "mean reversion volatilities"
%               rho: p-dimensional vector representing noise correlations
%                 t: T-dimensional vector of the time points for which
%                    a price is to be simulated
%
%    OUTPUT       y: T x n matrix of the simulated log prices
%            lambda: T x p matrix of the simulated eigenvalues
%
% See also SIM_WASC_2D.
%
% created by Benedikt Rudolph
% DATE: 12-Aug-2012
  
  n = size(A,1); % dimension of the problem
  p = size(A,2); % number of driving eigenvectors / noises
  T = length(t); % number of time grid points
  dt = reshape(diff(t),T-1,1); % time increments
  
  mu = reshape(mu, 1, n);
  kappa = reshape(kappa, 1, p);
  theta = reshape(theta, 1, p);
  sigma = reshape(sigma, 1, p);
  rho = reshape(rho, 1, p);
  
  y = zeros(T, n); % pre-allocation
  y(1,:) = y_0;
  
  lambda = zeros(T, p); % pre-allocation
  lambda(1,:) = lambda_0;
  
  dB = normrnd(0, 1, T-1, p); % noise driving the eigenvalues
  dW = normrnd(0, 1, T-1, p); % noise driving the log returns
  % correlate noises by rho
  dW = repmat(rho,T-1,1).*dB + repmat(sqrt(1-rho.^2),T-1,1).*dW;
  
%   for k=2:T
%     y(k,:) = y(k-1,:) + (mu - 0.5 * lambda(k-1,:) * A'.^2).*dt(k-1) ...
%                   + ( sqrt(lambda(k-1,:)) .* dW(k-1,:) * A' ) *sqrt(dt(k-1));
%     lambda(k,:) = lambda(k-1,:) + kappa.*(theta-lambda(k-1,:))*dt(k-1) ...
%                         + sigma.*sqrt(lambda(k-1,:)).*dB(k-1,:)*sqrt(dt(k-1));
%     lambda(k,:) = abs( lambda(k,:) );
%   end

  for k=2:T
    lambda(k,:) = lambda(k-1,:) + kappa.*(theta-lambda(k-1,:))*dt(k-1) ...
                        + sigma.*sqrt(lambda(k-1,:)).*dB(k-1,:)*sqrt(dt(k-1));
    lambda(k,:) = max( lambda(k,:),0 );%abs( lambda(k,:) );                
    y(k,:) = y(k-1,:) + (mu - 0.5 * lambda(k,:) * A'.^2).*dt(k-1) ...
                  + ( sqrt(lambda(k,:)) .* dW(k-1,:) * A' ) *sqrt(dt(k-1));               
  end



end







