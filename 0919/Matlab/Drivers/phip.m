function [w, stats] = phip(t, A, u, tol, symm, m)
% PHIP - Evaluates a linear combinaton of the phi functions
%        evaluated at tA acting on vectors from u, that is 
%
%        w = phi_0(tA) u(:, 1) + t phi_1(tA) u(:, 2) + 
%            t^2/2 phi_2(tA) u(:, 3) + ...  
%
%        The evaluation expresses eveything in terms of the highest
%        order phi function and evaluates the action of this on a
%        vector using a Krylov technique and then computes w using
%        the recurrence relation.
%
%        The size of the Krylov subspace is fixed during the
%        integration. The Krylov subspace is computed using Arnoldi if
%        A is non-symmetric and Lancozs if A is symmetric.
%
% PARAMETERS:
%   t    - constant value represent time.
%   A    - the matrix argument of the phi functions.
%   u    - the matrix with columns representing the vectors to be
%          multiplied by the phi functions.
%   tol  - the convergence tolarance required. 
%   symm - boolean true if matrix is symmetric.
%   m    - size of Krylov subspace.
%
% RETURNS:
%   w        - the linear combination of the phi functions
%              evaluated at tA acting on the vectors from u.
%   stats(1) - number of substeps
%   stats(2) - number of rejected steps
%   stats(3) - number of Krylov steps
%   stats(4) - number of matrix exponentials

persistent V int

% n is the size of the original problem
% p is the number of phi functions
[n, p] = size(u);
Aisahandle = isa(A, 'function_handle');
if ~Aisahandle
  nnze = nnz(A);
else
  nnze = 10*n; % wild guess
end;

% Add extra column of zeros if p=1
if p == 1
  p = 2;
  u = [u, zeros(size(u))];
end

% Check inputs
if nargin < 6
  m = min(n, 30);
  if nargin < 5
    if Aisahandle
      symm = false;
    elseif max(full(max(A-A'))) < 100*eps
      symm = true;
    else
      symm = false;
    end
    if nargin < 4
      tol = 1.0e-7;
      if nargin < 3
        error('phipm:NotEnoughInputs',...
              'Not enough input arguments.');
      end  
    end
  end
end

% Preallocate matrices
if isempty(V) && isempty(int)
  V = zeros(n, m+1); 
  int = zeros(n, p);
elseif numel(V) ~= n*(m+1)
  V = zeros(n, m+1); 
  int = zeros(n, p);  
elseif numel(int) ~= n*p
  int = zeros(n, p);  
end

% Initializing the variables
step = 0; 
krystep = 0;
ireject = 0;
reject = 0;
exps = 0;
happy = 0;
sgn = sign(t); 
tnow = 0;
tout = abs(t);
j = 0;

% Infinity norms
Anorm = norm(A, 'inf');
vnorm = norm(u(:, 1), 'inf');

% Compute and initial starting approximation for the timestep
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
tau = min(tout, (10/Anorm)*((fact*tol)/(4*vnorm*Anorm))^(1/m));

% Setting the safety factors and tolerance requirements
gamma = 0.8; 
delta = 1.2; 

% Used for toeplitz trick
cidx = (0:p-1)';
ridx = p:-1:1;
idx = cidx(:,ones(p,1))+ridx(ones(p,1),:);

% Initial condition
w = u(:, 1);
oldtau = NaN; omega = NaN; orderold = true;

% Iterate until we reach the final time t
while tnow < tout
  
  % Compute necessary starting information
  if j == 0
    
    % Initialize the matrices V and H
    H = zeros(m+p, m+p);
    x = [zeros(1, p-1), cumprod([1, tnow./(1:p-1)])];
    up = u*x(idx);  % Code inspired by toeplitz.m
        
    % Compute the update factors
    int(:, 1) =  w;
    for i = 1:p-1

      % Determine if matrix free
      if ~Aisahandle
        int(:, i+1) = A*int(:, i)+up(:, i+1);
      else
        int(:, i+1) = A(int(:, i))+up(:, i+1);
      end;
    
    end 

    % Normalize initial vector
    beta = norm(int(:, end));
    if beta == 0
    
      % Multiplying with a zero vector, hence result is zero
      % Finish all in one step
      reject = reject+ireject;
      step = step+1;
      tau = tout-tnow;
      w = w+int(:, 2:p-1)*cumprod(tau*1./(1: p-2)');
      break;
    
    end;
  
    % The first Krylov basis vector
    V(:, 1) = int(:, end)./beta;
  
  end;

  % Check if matrix is symmetric 
  if symm
  
    % Symmetric use the Lanczos process
    while j < m

      % Determine if matrix free
      j = j+1;
      if ~Aisahandle
        vv = A*V(:, j);
      else
        vv = A(V(:, j));
      end;
      H(j, j) = V(:, j)'*vv;     
      if j == 1,
        vv = vv-H(j, j)*V(:, j); 
      else
        vv = vv-H(j-1, j)*V(:, j-1)-H(j, j)*V(:, j); 
      end       
      krystep = krystep+1;
      s = norm(vv);
      
      % Happy breakdown
      if s < tol
        happy = 1;
        tau = tout-tnow;
        break;
      end;
      H(j+1, j) = s;
      H(j, j+1) = s;
      V(:, j+1) = vv./s;
    
    end;

    % Keep a record of H
    H2 = H;
    H(m, m+1) = 0; 
  
  % Matrix is not symmetric  
  else
  
    % Not symmetric use the Arnoldi process
    while j < m

      % Determine if matrix free
      j = j+1;
      if ~Aisahandle
        vv = A*V(:, j);
      else
        vv = A(V(:, j));
      end;
      for i = 1:j
        H(i, j) = V(:, i)'*vv;
        vv = vv-H(i, j)*V(:, i);
      end
      krystep = krystep+1;
      s = norm(vv);
      
      % Happy breakdown
      if s < tol
        happy = 1;
        tau = tout-tnow;
        break;
      end;
      H(j+1, j) = s;
      V(:, j+1) = vv./s;
    
    end 
  
    % Keep a record of H
    H2 = H;
  
  end

  % We use the vector e1 in the computations
  H(1, m+1) = 1; 
  
  % Construct the augmented matrix
  for i = 1:p-1
    H(m+i, m+i+1) = 1;
  end
  h = H(m+1, m); 
  H(m+1, m) = 0;
  
  % Compute the exponential of the augmented matrix
  %[F, hnorm] = expmhigham(sgn*tau*H);
  F = expm(sgn*tau*H);
  exps = exps+1;
  
  % Local truncation error estimation
  err = abs(beta*h*F(m, m+p));
  
  % Error per unit step
  oldomega = omega;
  omega = tout*err/(tau*tol);

  % Estimate the order
  if tau ~= oldtau && ireject >= 1
    order = max(1, log(omega/oldomega)/log(tau/oldtau));
    orderold = false;
  elseif orderold || ireject == 0
    orderold = true;
    order = j/4;
  else
    orderold = true;
  end;
  
  % Update tau
  oldtau = tau;
  taunew = tau*(omega/gamma)^(-1/order);

  % Check error against target
  if omega <= delta   

    % Yep, got the required tolerance; update
    reject = reject+ireject;
    step = step+1;
    up = w+int(:, 2:end-1)*cumprod(tau*1./(1: p-2)');
    % up = w+sum_{k=1}^{p-2} tau^k / k! * int(:,k+1)
    
    % Using the corrected quantity
    F(m+1, m+p-1) = h*F(m, m+p);
    w = beta*V(:, 1:m+1)*F(1:m+1, m+p-1)+up;
    
    % Update t
    tnow = tnow+tau;
    j = 0;
    ireject = 0;
    
  else
    
    % Nope, try again
    H = H2;
    ireject = ireject+1;
  
  end;
  
  % Safety factors for tau
  tau = min(tout-tnow, max(tau/5, min(5*tau, taunew)));
    
end
stats = [step, reject, krystep, exps];
