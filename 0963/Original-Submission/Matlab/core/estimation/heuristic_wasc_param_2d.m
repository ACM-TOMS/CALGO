function [mu, Sigma_0, M, Q, rho, beta] = heuristic_wasc_param_2d(y_t, dt)
%HEURISTIC_WASC_PARAM Obtains heuristic parameter estimates for a WASC model.
%
%  [mu, Sigma_0, M, Q, rho, beta] = heuristic_pcsv_param(y_t, dt)
%    calculates parameter estimates for a WASC model with the time series
%    y_t using a heuristic approach. The results are supposed to be used
%    as starting points for an iterative optimization algorithm.
%
%    INPUT y_t: A Txn matrix representing the log prices
%           dt: The time difference in the observed series y_t
%
% created by Benedikt Rudolph
% DATE: 17-Oct-2012

  [T, n] = size(y_t);
  r = diff(y_t);
  % Obtain mu and Sigma_0 == Sigma_infty
  Sigma_0 = cov(r)/dt;
  mu = mean(r)'/dt + (0.5*diag(Sigma_0));
  % Calculate a rolling time series of covariance matrices as a basis
  % for the calculation of Q, M and beta
  roll_period = 25;
  S = zeros(n, n, T-roll_period);
  for k=1:(T-roll_period)
    S(:,:,k) = cov(diff(y_t(k:(k+roll_period-1),:)))/dt;
  end
  % Calculate moments of dS and S
  dS = diff(S, 1, 3); % first-order difference along third dimension
  K = size(dS,3);
  E = @(a) a'*mean(dS,3)*a;
  V = @(a) mean( (arrayfun(@(k) a'*dS(:,:,k)*a, 1:K) - E(a)).^2 );
  Cov = @(a,g) mean( ...
                (arrayfun(@(k) a'*dS(:,:,k)*a, 1:K) - E(a)) ...
            .*  (arrayfun(@(k) g'*dS(:,:,k)*g, 1:K) - E(g)) );
  KS = size(S,3);
  ES = @(a) a'*mean(S,3)*a;
  VS = @(a) mean( (arrayfun(@(k) a'*S(:,:,k)*a, 1:KS) - ES(a)).^2 );
  K = @(a) 2*ES(a)^2 / VS(a);
  a1 = [1;0];
  a2 = [0;1];
  a12 = [1;1];
  % Obtain beta
  beta = K(a12);
  beta = max(beta, 9);
  % Obtain Q
  Q = zeros(n);
  Q(1,1) = sqrt( V(a1) / (4*Sigma_0(1,1)*dt) );
  Q(2,2) = sqrt( V(a2) / (4*Sigma_0(2,2)*dt) );
  Q(1,2) = -Cov(a1,a2) / (4*Sigma_0(1,2)*Q(1,1)*dt);
  % Obtain M as solution to A*M = b
  A = zeros(3);
  A(1,1) = 2*Sigma_0(1,1);
  A(1,2) = 2*Sigma_0(1,2);
  A(2,2) = 2*Sigma_0(1,2);
  A(2,3) = 2*Sigma_0(2,2);
  A(3,1) = 2*(Sigma_0(1,1)+Sigma_0(1,2));
  A(3,2) = 2*(Sigma_0(1,1)+2*Sigma_0(1,2)+Sigma_0(2,2));
  A(3,3) = 2*(Sigma_0(2,2)+Sigma_0(1,2));
  b = zeros(3,1);
  b(1) = E(a1)/dt - beta*a1'*(Q'*Q)*a1;
  b(2) = E(a2)/dt - beta*a2'*(Q'*Q)*a2;
  b(3) = E(a12)/dt - beta*a12'*(Q'*Q)*a12;
  m = A \ b;
  M = [m(1), m(2); m(2), m(3)];
  % Obtain rho
  period_length = 50;
  number_of_periods = floor((T-1)/period_length)-1;
  dZ = zeros(number_of_periods*n, 1);
  dW = zeros(number_of_periods*n, n);
  S_last = Sigma_0;
  for k=1:number_of_periods
    idx = (k*period_length):((k+1)*period_length);
    period_returns = r(idx, :);
    S = cov(period_returns);
    dS = (S - S_last) / period_length;
    dY = mean(period_returns);
    sqrt_S = chol(S);
    dZ((n*(k-1)+1):(n*k)) = sqrt_S \ (dY' - (mu - diag(S))*dt);
    dW((n*(k-1)+1):(n*k), :) = (2*sqrt_S*Q) \ (dS-(beta*(Q'*Q)+2*M*S)*dt);
  end
  rho = (dW \ dZ);
end
