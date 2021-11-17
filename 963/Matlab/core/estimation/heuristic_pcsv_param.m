function [mu, A, lambda_0, kappa, theta, sigma, rho] ...
          = heuristic_pcsv_param(y_t, dt, p)%, mu1, A1, lambda_01)
%function [mu, A, lambda_0, kappa, theta, sigma, rho] ...
%          = heuristic_pcsv_param(y_t, dt, p)
%HEURISTIC_PCSV_PARAM Obtains heuristic parameter estimates for a PCSV model.
%
%  [mu, A, lambda_0, kappa, theta, sigma, rho]
%    = heuristic_pcsv_param(y_t, dt) calculates parameter estimates for a PCSV
%    model with the time series y_t using a heuristic approach. The results
%    are supposed to be used as starting points for an iterative optimization
%    algorithm.
%
%    INPUT y_t: A Txn matrix representing the log prices
%           dt: The time difference in the observed series y_t
%            p: Number of eigenvectors used for this model
%
% created by Benedikt Rudolph
% DATE: 16-Aug-2012

  [T, n] = size(y_t);
  
  % obtain A and lambda_0
  [A, Lambda, V] = svd(cov(diff(y_t))/dt);
  
  lambda_0 = diag(Lambda);
  
  % truncate to desired number of eigenvectors/eigenvalues
  lambda_0 = lambda_0(1:p)';
  A = A(:,1:p);
  
  % obtain mu
  mu = (mean(diff(y_t))/dt + (0.5*(A.^2)*lambda_0')')';
  
  % obtain a time series l of eigenvalues using rolling estimates
  roll_period = 25;%8 new
  l = zeros(T-roll_period, p);
  %_dY = zeros(T-roll_period, p);
  r_t = diff(y_t);
  %dB = zeros(T-roll_period-1, p);
  %dW = zeros(T-roll_period-1, p);
  for k=1:(T-roll_period)
    S = cov(r_t(k:(k+roll_period-1),:))/dt;% cov(y_t(k:(k+roll_period-1),:))/dt;
    %[A_k, Lambda_k, V_k] = svd(S);
    lambda_k = diag(A'*S*A); % diag(Lambda_k); %(A'*r_t(k+roll_period-1,:)').^2; %
    l(k,:) = lambda_k(1:p); % /dt;
    %_dY(k,:) = mean(diff(y_t(k:(k+roll_period-1),:)));
    %dW(k,:) = (A*diag(sqrt(l(k,:))))*sqrt(dt) ...
    %         \ (r_t(k,:)' - (mu - 0.5*A.^2*l(k,:)')*dt);
  end
  
  % Obtain kappa,theta and sigma using OLS on the eigenvalue estimates l
  kappa = zeros(1,p);
  theta = zeros(1,p);
  sigma = zeros(1,p);
  feller_condition = @(kappa, theta, sigma) 2*kappa.*theta > sigma.^2;
  for k=1:p
    lambda = l(:,k);
    Y = diff(lambda)./sqrt(lambda(1:(end-1)));
    X = [dt./sqrt(lambda(1:(end-1))), dt*sqrt(lambda(1:(end-1)))];
    b = X \ Y; %b=lscov(X,Y)
    r = Y - X*b;
    s = (r'*r) / (size(X,1)-2);
    %dB(:,k) = r/sqrt(s);
    sigma(k) = sqrt(s/dt);
    %theta(k) = -b(1) / b(2);
    theta(k) = lambda_0(k);
    kappa(k) = -b(2);
    % if kappa, theta and sigma do not satisfy the Feller condition
    % shift the parameters equally so that they do satisfy it
    if ~feller_condition(kappa(k), theta(k), sigma(k))
      alpha = - ( sqrt(2*kappa(k).*theta(k)) - sigma(k) ) ./ ...
                ( sqrt(2*kappa(k).*theta(k)) + sigma(k) ) * 1.01;
      kappa(k) = (1+alpha) * kappa(k);
      theta(k) = (1+alpha) * theta(k);
      sigma(k) = (1-alpha) * sigma(k);
    end
  end
  
   % Obtain rho
  period_length = 30;%12;
  number_of_periods = floor((T-1)/period_length)-1;
  dB = zeros(number_of_periods, p);
  dW = zeros(number_of_periods, p);
  lambda_last = lambda_0;
  for k=1:number_of_periods
    idx = (k*period_length):((k+1)*period_length);
    period_returns = r_t(idx, :);
    S = cov(period_returns);
    %[A_k, Lambda_k, V_k] = svd(S);
    lambda_k = diag(A'*S*A); %diag(Lambda_k);
    lambda_k = lambda_k(1:p);
    dlambda = (lambda_k'-lambda_last) / period_length;
    dY = mean(period_returns);
    dB(k,:) = ( (dlambda - kappa.*(theta-lambda_k')*dt) ./ ...
                                  (sigma .* sqrt(lambda_k')) ) / sqrt(dt);
    dW(k, :) = ( (A*diag(sqrt(lambda_k)))*sqrt(dt) ...
              \ (dY' - (mu - 0.5*A.^2*lambda_k)*dt) )' / sqrt(dt);
  end
  rho = diag(corr(dW,dB))';
end
