function phi = cf_wasc_v(mu, Sigma_0, M, Q, rho, beta ...
                              , omega, y_t, tau)
% CF_WASC_V Calculates characteristic function values for a WASC model.
%
%  phi = cf_wasc_v(mu, Sigma_0, M, Q, rho, beta, omega, y_t, tau)
%    calculates the vectorized conditional characteristic function
%    of continuous returns of a p-dimensional WASC model as described
%    in da Fonseca et al. 2007 (all references to this paper).
%    Note that the cf in eq. (21) does not depend on y_t, but it is included
%    in the argument list for compatibility reasons (an arbitrary value
%    can be passed).
%    Assuming that there are T observations (=time points),
%    N grid points for the cf evaluation,
%    p the dimension of the problem,
%    the parameters accepted will be:
%
%    INPUT      mu: a p-dimensional mean vector mu
%          Sigma_0: a p x p covariance matrix
%                M: a p x p matrix
%                Q: an orthogonal p x p matrix
%              rho: a p-dimensional correlation vector
%             beta: a real value
%            omega: a p x N matrix containing the evaluation grid
%              y_t: a T x p matrix containing the observations (only
%                   the length is used to determin T)
%              tau: a real valued time difference
%
%    OUTPUT    phi: a (T-1) x N matrix epresenting the cf value
%                   for time (row) and grid point (column).
%
% See also CF_PCSV_V, CF_GBM_V.
%
% created by Benedikt Rudolph
% DATE: 20-Aug-2012
  
  % p: number of assets
  % N: number of grid points
  [p, N] = size(omega);
  T = size(y_t,1);
  t = cumsum(repmat(tau,T,1))-tau;
  
  phi = zeros(T-1, N);% pre-allocate phi
  
  A_large = compute_a_large(M, Q, rho, omega, tau);
  B_large = compute_b_large(M, Q, t(1:end-1)); % note call with t for tau! eq. (21)
  trace_M = trace(M);
  t_trace_M = trace_M*t(1:end-1);
  
  i_1_to_p = 1:p;
  parfor k=1:N
    A_21 = A_large(p+(1:p),2*p*(k-1)+(1:p)); % eq. (14)
    A_22 = A_large(p+(1:p),2*p*(k-1)+p+(1:p)); % eq. (14)
    A = A_22 \ A_21; % eq. (13)
    omega_k = omega(:,k);
    c = -0.5*beta*(real(log(det(A_22)))+tau*trace_M ...
      + tau*1i*sum(sum((omega_k*rho').*Q'))) ...
      + tau*1i*mu'*omega_k;  % eq. (15)
    Delta = -1i*A;  % eq. (21)
    phi_latest = 1;
    for j=1:(T-1)
      i_2pj_1 = 2*p*(j-1);
      B_11 = B_large(i_2pj_1+i_1_to_p,i_1_to_p);
      B_12 = B_large(i_2pj_1+i_1_to_p,p+i_1_to_p);
      B_21 = B_large(i_2pj_1+p+i_1_to_p,i_1_to_p);
      B_22 = B_large(i_2pj_1+p+i_1_to_p,p+i_1_to_p);
      to_inv = 1i*Delta*B_12+B_22;
      if ~( rcond(to_inv)<sqrt(eps) )
        B = to_inv \ (1i*Delta*B_11+B_21); % eq. (18)
        trace_B_Sigma_0 = sum(sum(B.*Sigma_0'));
        trace_log_D_B = real(log(det(to_inv)));
        % note call with t for tau (!) eq. (21)
        C = -0.5*beta*(trace_log_D_B+t_trace_M(j)); % eq. (18)
        phi_latest = c+trace_B_Sigma_0+C; % eq. (21)/(17)
      end
      phi(j,k) = phi_latest;
    end
  end
  phi = exp(phi);
end

function A_large = compute_a_large(M, Q, rho, omega, tau)
  % Calculates the 2p x 2p matrix described in eq. (14).
  % That matrix is referred to as "A_large" here in the code.
  % The vectorized version of A_large is a 2pT x 2pN matrix. 
  %
  % A_large = compute_a_large(M, Q, rho, omega, tau)
  
  p = size(omega, 1); % number of assets == dimension of the problem
  N = size(omega, 2);
  
  % pre-calulate matrix elements of the matrix to exponentiated in eq. (14)
  
  omega_2 = zeros(p, p*N); % pre-allocate
  % compute matrix omega_k * omega_k' for each grid point
  for k=1:N
    omega_2(:,(k-1)*p+1:k*p) = omega(:,k) * omega(:,k)';
  end
  
  omega_flat = reshape(omega, 1, p*N); % flatten omega to a single row vector
  
  ul = repmat(M, 1, N) + 1i*Q' * rho * omega_flat; % p x p*N matrix
  ll = 0.5*(-omega_2 - 1i*repmat(eye(p),1,N).*repmat(omega_flat,p,1)); % p x p*N matrix
  ur = -2*(Q'*Q); % p x p matrix
  lr = -(repmat(M',1,N) + 1i*kron(omega,rho'*Q)); % p x p*N matrix
    
  % distribute the four parts of the matrix to be exponentiated
  % on a 2*p x 2*p*N matrix
  a = zeros(2*p, 2*p*N);
  a(repmat(logical([ones(p),zeros(p);zeros(p,2*p)]),1,N)) = ul;
  a(repmat(logical([zeros(p),ones(p);zeros(p,2*p)]),1,N))=repmat(ur,1,N);
  a(repmat(logical([zeros(p,2*p);ones(p),zeros(p)]),1,N)) = ll;
  a(repmat(logical([zeros(p,2*p);zeros(p),ones(p)]),1,N)) = lr;
  
  % exponentiate
  p2 = 2*p;
  A_large_cell = cell(1,N);
  parfor j=1:N
    A_large_cell{j} = expm(tau*a(:,(j-1)*p2+1:j*p2));
  end
  
  A_large = [A_large_cell{:}]; % 2*p x 2*p*N matrix
end

function B_large = compute_b_large(M, Q, tau)
  % Calculates the 2p x 2p matrix described in eq. (18) in a vectorized way.
  % That matrix is referred to as "B_large" here in the code.
  % The vectorized version of B_large is a 2p(T-1) x 2p matrix. 
  %
  % B_large = compute_b_large(M, Q, tau)
  
  p = size(M, 1); % number of assets == dimension of the problem
  T = length(tau);
  
  B_log = [M,-2*(Q'*Q);zeros(p),-M']; % matrix of which exp is to be computed
  
  p2 = 2*p;
  B_large = zeros(p2*(T-1), p2);
  
  for j=1:T
    B_large((j-1)*p2+1:j*p2,:) = expm(tau(j)*B_log);
  end
  
end

