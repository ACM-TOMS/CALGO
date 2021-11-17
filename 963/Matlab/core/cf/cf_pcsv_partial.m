function phi = cf_pcsv_partial(mu, A, lambda_0, kappa, theta, sigma ...
                              , rho, omega, y_t, tau)
%CF_PCSV_PARTIAL Calculates characteristic function values for a PCSV partial model.
%
%  phi = cf_pcsv_v(mu, A, lambda_0, kappa, theta, sigma, rho, omega, y_t, tau)
%    Calculates the conditional characteristic function of continuous returns
%    of an n-dimensional PCSV model in a vectorized way.
%    Assuming that n is the dimension of the problem, p the number of
%    eigenvalues/eigenvectors, the parameters expected are
%
%    INPUT       mu: n-dimensional mean vector
%                 A: n x p matrix of eigenvectors
%          lambda_0: p-dimensional vector of initial eigenvalues
%             kappa: "mean reversion speed"
%             theta: "mean reversion means"
%             sigma: "mean reversion volatilities"
%               rho: noise correlation
%             omega: n x N matrix representing the evaluation grid
%               y_t: T x n matrix representing the time grid
%               tau: real valued time difference
%
%    OUTPUT     phi: a (T-1) x N matrix representing the cf value
%                    for time (row) and grid point (column).
%
% See also CF_PCSV_V, CF_GBM_V.
%
% created by Benedikt Rudolph
% DATE: 31-Jan-2013
  
  % n: number of assets
  % p: number of eigenvectors, only the 1st is stochastic
  [n, p] = size(A);
  
  N = size(omega, 2); % number of (n-dimensional) evaluation grid points
  T = size(y_t,1)-1; % number of time grid points
  t = cumsum(repmat(tau,T,1))-tau;
  
  c = omega'*A; % N x p
  b = -0.5*omega'*(A.^2); % N x p
  r = omega'*mu/p; % N x 1
  r = repmat(r, 1, p); % N x p
  
  ph = 1; % suppose we do not need it
  
  % reshape lambda_0 to required form
  lambda_0 = reshape(lambda_0, 1, p);
  
  lambda_0 = repmat(lambda_0, N, 1); % N x p
  kappa = repmat(kappa, N, 1); % N x 1
  theta = repmat(theta, N, 1); % N x 1
  sigma = repmat(sigma, N, 1); % N x 1
  rho = repmat(rho, N, 1); % N x 1
  
  % eigenvalues 2 to p
  
  phi_lognormal = prod( (exp(...
      (r(:,2:p) + b(:,2:p).*lambda_0(:,2:p))*tau*ph ...
      - 0.5 * (c(:,2:p).^2.*lambda_0(:,2:p))*tau*ph^2 ...
  )).', 1);
  
  % first eigenvalue
  
  c = c(:,1);
  b = b(:,1);
  r = r(:,1);
  
  d = sqrt((kappa-rho.*c.*sigma.*ph*1i).^2-sigma.^2.*(2*b.*ph*1i-c.^2.*ph.^2));
  g = (kappa - rho.*c.*sigma.*ph*1i + d) ./ (kappa - rho.*c.*sigma.*ph*1i - d);
  
  C = r*tau.*ph*1i + kappa.*theta./(sigma.^2) .* ...
      ((kappa-rho.*c.*sigma.*ph*1i+d)*tau - 2*log((1-g.*exp(d*tau))./(1-g)));
  D = (kappa - rho.*c.*sigma.*ph*1i + d) ./ (c.^2.*sigma.^2) .* ...
      (1 - exp(d*tau)) ./ (1-g.*exp(d*tau));
  
  Delta = -1i*D;
  
  phi = zeros(T,N);
  
  parfor k=1:length(t)
    tau = t(k);
 
    A = -2*kappa.*theta./(sigma.^2) ...
        .* log(1-1i*Delta.*c.^2.*sigma.^2./(2*kappa).*(1-exp(-kappa*tau))); 
    B = 1i*Delta.*exp(-kappa*tau) ...
        ./ (1 - c.^2.*sigma.^2*1i.*Delta./(2*kappa) .* (1-exp(-kappa*tau))); 

    phi(k,:) = phi_lognormal .* exp(C + A + B.*c.^2.*lambda_0(:,1)).';
  end
  
  % omega==0 yields NaN
  idx = find(sum(omega==0)==n); % omega(:,idx)==zeros(n,1)
  if idx
    phi(:,idx) = 1; % Phi(0) == 1 for all cf
  end
end
