function [y, lambda] = sim_pcsv_partial(y_0, mu, A, lambda_0, kappa ...
                                       , theta, sigma, rho, t)
%SIM_PCSV_PARTIAL Simulates a PCSV partial model for the given time points in t.
%
%  [y, lambda] = sim_pcsv_partial(y_0, mu, A, lambda_0, kappa, theta, sigma, rho, t)
%    Performs an Euler simulation of the PCSV partial model specified by the
%    given parameters for the time points in t.
%    Assuming that n is the dimension of the problem, p the number of
%    eigenvalues/eigenvectors, the parameters expected are
%
%    INPUT      y_0: n-dimensional vector with initial log prices
%                mu: n-dimensional mean vector
%                 A: n x p matrix of eigenvectors
%          lambda_0: p-dimensional vector of initial/constant eigenvalues
%             kappa: scalar "mean reversion speed"
%             theta: scalar "mean reversion mean"
%             sigma: scalar "mean reversion volatility"
%               rho: scalar noise correlation
%                 t: T-dimensional vector of the time points for which
%                    a price is to be simulated
%
%    OUTPUT       y: T x n matrix of the simulated log prices
%            lambda: T-dimensional vector of the simulated eigenvalues
%
% See also SIM_WASC_2D, SIM_PCSV.
%
% created by Benedikt Rudolph
% DATE: 04-Feb-2013
  
  n = size(A,1); % dimension of the problem
  p = size(A,2); % number of driving eigenvectors / noises, just one stochastic
  T = length(t); % number of time grid points
  dt = reshape(diff(t),T-1,1); % time increments
  
  mu = reshape(mu, 1, n);
  
  y = zeros(T, n); % pre-allocation
  y(1,:) = y_0;
  
  lambda = zeros(T, 1); % pre-allocation
  lambda(1) = lambda_0(1);
  
  dB = normrnd(0, 1, T-1, 1); % noise driving the eigenvalue
  dW = normrnd(0, 1, T-1, p); % noise driving the log returns
  % correlate noises by rho
  dW(:,1) = repmat(rho,T-1,1).*dB + repmat(sqrt(1-rho.^2),T-1,1).*dW(:,1);
  
  for k=2:T
    l_k = [lambda(k-1), lambda_0(2:end)];
    y(k,:) = y(k-1,:) + (mu - 0.5 * l_k * A'.^2).*dt(k-1) ...
                  + ( sqrt(l_k) .* dW(k-1,:) * A' ) *sqrt(dt(k-1));
    lambda(k) = lambda(k-1) + kappa.*(theta-lambda(k-1))*dt(k-1) ...
                        + sigma.*sqrt(lambda(k-1)).*dB(k-1)*sqrt(dt(k-1));
    lambda(k,:) = abs( lambda(k,:) );
  end
end
