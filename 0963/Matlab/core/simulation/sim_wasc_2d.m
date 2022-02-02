function [y, vol] = sim_wasc_2d(y0, mu, Sigma_0, M, Q, rho, beta, t)
%SIM_WASC_2D Simulates a 2d WASC model for the given time points in t.
%
%  [y, vol] = sim_wasc_2d(y0, mu, Sigma_0, M, Q, rho, beta, t)
%    Performs an Euler simulation of the WASC model specified by the given
%    parameters for the time points in t.
%    Assuming that n is the dimension of the problem, the parameters expected are
%
%    INPUT      y_0: n-dimensional vector with initial log prices
%                mu: n-dimensional mean vector mu
%           Sigma_0: n x n covariance matrix
%                 M: n x n matrix
%                 Q: orthogonal n x n matrix
%               rho: n-dimensional correlation vector
%              beta: a natural number
%                 t: T-dimensional vector of the time points for which
%                    a price is to be simulated
%
%    OUTPUT       y: T x n matrix of the simulated log prices
%               vol: T x 3 matrix of the simulated variances and correlation
%
% See also SIM_PCSV.
%
% created by Benedikt Rudolph
% DATE: 12-Aug-2012
  
  mu = reshape(mu, 1, 2);
  T = length(t); % number of time grid points
  dt = reshape(diff(t),T-1,1); % time increments
  
  y = zeros(T, 2); % pre-allocation
  y(1,:) = y0;
  
  vol = zeros(T, 3);
  vol(1,:) = [ sqrt(Sigma_0(1,1)), sqrt(Sigma_0(2,2)), ...
    Sigma_0(1,2)/sqrt(Sigma_0(1,1)*Sigma_0(2,2)) ];
  
  % choose values for X_0 using eigenvalues and
  % eigenvectors of Sigma_0
  [V, Lambda] = eig(Sigma_0);
  X = repmat(sqrt(diag(Lambda)),1,2).*V';
  X = [X; zeros(beta-2, 2)]; % add zeros in the remaining part
  
  % init matrices which are updated in every step
  Sigma = Sigma_0;
  sqrtm_Sigma = sqrtm(Sigma);
  sigma_tilde = inv(sqrtm_Sigma);
  
  for k=2:T
    
    % simulate noise
    dN = sqrt(dt(k-1))*normrnd(0,1,beta,2);
    dB = sqrt(dt(k-1))*normrnd(0,1,2,1);
    dW = sigma_tilde*X'*dN;
    
    % Euler approximation of y
    y(k,:) = y(k-1,:) + (mu - 0.5*diag(Sigma)')*dt(k-1) ...
                      + (sqrtm_Sigma*(dW*rho + sqrt(1-rho'*rho)*dB))';
    
    X = X + X*M*dt(k-1) + dN*Q; % Euler approximation of X
    % updates of Sigma
    Sigma = X'*X;
    sqrtm_Sigma = sqrtm(Sigma);
    sigma_tilde = inv(sqrtm_Sigma);

    % store volatilites and correlation
    v1 = sqrt(Sigma(1,1));
    v2 = sqrt(Sigma(2,2));
    vol(k,:) = [ v1, v2, Sigma(1,2)/(v1*v2) ];
  end
end
