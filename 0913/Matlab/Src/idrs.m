function [x,flag,relres,iter,resvec,replacements]=idrs(A,b,s,tol,maxit,M1,M2,x0,options );
%IDRS Induced Dimension Reduction method
%   X = IDRS(A,B) solves the system of linear equations A*X=B for X.  
%   The N-by-N coefficient matrix A must be square and the right-hand
%   side column vector B must have length N. A can also be a struct. The field 
%   A.name should contain the name of a function to perform multiplications with A
%   Other fields can be used to pass parameters to this function.    
%
%   X = IDRS(A,B,S) specifies the dimension of the 'shadow space'. If S = [], then
%   IDRS uses the default S = 4. Normally, a higher S gives faster convergence, 
%   but also makes the method more expensive.
%
%   X = IDRS(A,B,S,TOL) specifies the tolerance of the method.  If TOL is []
%   then IDR uses the default, 1e-8.
%
%   X = IDRS(A,B,S,TOL,MAXIT) specifies the maximum number of iterations.  If
%   MAXIT is [] then IDRS uses the default, min(2*N,1000).
%
%   X = IDRS(A,B,S,TOL,MAXIT,M1) use preconditioner M1. If M1 is [] then no 
%   preconditioner is applied. 
%   IDRS(A,B,S,TOL,MAXIT,M1,M2) uses a factored preconditioner M = M1 M2. 
%   M1 and M2 can be structures. In that case the field M1.name and M2.name 
%   should contain function names for M1 and M2.
%
%   X = IDRS(A,B,S,TOL,MAXIT,M1,M2,X0) specifies the initial guess.  If X0 is []
%   then IDR uses the default, an all zero vector.
%
%   X = IDRS(A,B,S,TOL,MAXIT,M1,M2,X0,OPTIONS) specifies additional options.
%      OPTIONS must be a structure 
%      OPTIONS.SMOOTHING specifies if residual smoothing must be applied
%         OPTIONS.SMOOTHING = 0: No smoothing
%         OPTIONS.SMOOTHING = 1: Smoothing
%         Default: OPTIONS.SMOOTHING = 0;
%      OPTIONS.OMEGA determines the computation of OMEGA
%         If OPTIONS.OMEGA = 0: a standard minimum residual step is performed
%         If OPTIONS.OMEGA > 0: OMEGA is increased if
%         the cosine of the angle between Ar and r < OPTIONS.OMEGA
%         Default: OPTIONS.OMEGA = 0.7;
%      OPTIONS.P defines the 'shadow' space
%         Default: OPTIONS.P = ORTH(RANDN(N,S));
%      OPTIONS.REPLACE determines the residual replacement strategy
%         If |r| > 1E3 |b| TOL/EPS) (EPS is the machine precision)
%         the recursively computed residual is replaced by the true residual
%         once |r| < |b| (to reduce the effect of large intermediate residuals 
%         on the final accuracy)
%         Default: OPTIONS.REPLACE = 0; (No residual replacement)
%
%   [X,FLAG] = IDRS(A,B,S,TOL,MAXIT,M1,M2,X0,OPTIONS) 
%   also returns an information flag:
%       FLAG = 0: required tolerance satisfied
%       FLAG = 1: no convergence to the required tolerance within maximum 
%                 number of iterations
%       FLAG = 2: check RELRES, possible stagnation above required 
%                 tolerance level
%       FLAG = 3: one of the iteration parameters became zero, causing break down
%
%   [X,FLAG,RELRES] = IDRS(A,B,S,TOL,MAXIT,M1,M2,X0,OPTIONS) also 
%   returns the relative residual norm:
%          RELRES = ||B - AX||/||B||
%   
%   [X,FLAG,RELRES,ITER] = IDRS(A,B,S,TOL,MAXIT,M1,M2,X0,OPTIONS) also returns
%   the number of iterations.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = IDRS(A,B,S,TOL,MAXIT,M1,M2,X0,OPTIONS) also returns
%   a vector of the residual norms at each matrix-vector multiplication.
%
%   [X,FLAG,RELRES,ITER,RESVEC,REPLACEMENTS] = IDRS(A,B,S,TOL,MAXIT,M1,M2,X0,OPTIONS) 
%   also returns the number of residual replacements
%
%   Martin van Gijzen 
%   Version August 31, 2010
%   
%   This software is distributed under the
%   ACM Software Copyright and License Agreement. 
%

if ( nargout == 0 )
   help idrs;
   return
end

% Check for an acceptable number of input arguments
if nargin < 2
   error('Not enough input arguments.');
end

% Check matrix and right hand side vector inputs have appropriate sizes
funA = 0;
if isa(A,'struct') 
   funA = 1;
   if isfield(A,'name')
      function_A = A.name;
   else
      error('Use field A.name to specify function name for matrix-vector multiplication');
   end
   n = length(b);
else
   [m,n] = size(A);
   if (m ~= n)
      error('Matrix must be square.');
   end
   if ~isequal(size(b),[m,1])
      es = sprintf(['Right hand side must be a column vector of' ...
            ' length %d to match the coefficient matrix.'],m);
      error(es);
   end
end

% Assign default values to unspecified parameters
if nargin < 3 || isempty(s)
   s = 4;
end
if ( s > n )
   s = n;
end
if nargin < 4 || isempty(tol)
   tol = 1e-8;
end
if nargin < 5 || isempty(maxit)
   maxit = min(2*n,1000);
end

if nargin < 6 || isempty(M1)
   precL = 0;
elseif isa(M1,'struct') 
   if isfield(M1,'name')
      function_M1 = M1.name;
   else
      error('Use field M1.name to specify function name preconditioner');
   end
   precL = 1;
   funL = 1;
else
   if ~isequal(size(M1),[n,n])
      es = sprintf(['Preconditioner must be a matrix of' ...
            ' size %d times %d to match the problem size.'],n,n);
      error(es);
   end
   precL = 1;
   funL = 0;
end

if nargin < 7 || isempty(M2)
   precU = 0;
elseif isa(M2,'struct') 
   if isfield(M2,'name')
      function_M2 = M2.name;
   else
      error('Use field M2.name to specify function name preconditioner');
   end
   precU = 1;
   funU = 1;
else
   if ~isequal(size(M2),[n,n])
      es = sprintf(['Preconditioner must be a matrix of' ...
            ' size %d times %d to match the problem size.'],n,n);
      error(es);
   end
   precU = 1;
   funU = 0;
end

if nargin < 8 || isempty(x0)
   x0 = zeros(n,1);
else
   if ~isequal(size(x0),[n,1])
      es = sprintf(['Initial guess must be a column vector of' ...
            ' length %d to match the problem size.'],n);
      error(es);
   end
end

% Other parameters 
smoothing = 0;
angle = 0.7;
replacement = 0;
replacements = 0;
randn('state', 0);
P = randn(n,s);
P = orth(P);

if ( nargin > 8  ) 
% Residual smoothing:
   if isfield(options,'smoothing') 
      smoothing = options.smoothing > 0;
   end
% Computation of omega:
   if isfield(options,'omega')
      angle = options.omega;
   end 
% Alternative definition of P:
   if isfield(options,'P')
      P = options.P;
      if ~isequal(size(P),[n,s])
         es = sprintf(['P must be a matrix of' ...
         ' size %d times %d to match the problem size.'],n,s);
         error(es);
      end
   end 
   if isfield(options,'replace' )
      replacement = options.replace > 0;
   end 
end

if nargin > 9
   es = sprintf(['Too many input parameters']);
   error(es);
end

% END CHECKING INPUT PARAMETERS AND SETTING DEFAULTS

% Check for zero rhs:
if (norm(b) == 0)              % Solution is nulvector
   x = zeros(n,1);        
   iter = 0;                 
   resvec = 0;
   flag = 0;
   relres = 0;
   return
end
%
% Number close to machine precision:
mp = 1e3*eps;
%
% Initialize output paramater relres
relres = NaN;
%
% Compute initial residual:
x = x0;
normb = norm(b);
tolb = tol * normb;           % Relative tolerance
%
if funA
   r = b- feval( function_A, x, A);
else 
   r = b - A*x;
end
%
if smoothing
   x_s = x0;
   r_s = r;
end
%
normr = norm(r);
resvec=[normr];
trueres = 0;
%
if (normr <= tolb)           % Initial guess is a good enough solution
   iter = 0;                 
   flag = 0;
   relres = normr/normb;
   return
end
%
G = zeros(n,s); U = zeros(n,s); M = eye(s,s); 
om = 1;
%
% Main iteration loop, build G-spaces:
iter = 0;
while ( normr > tolb && iter < maxit )  
%
% New righ-hand size for small system:
   f = (r'*P)';
   for k = 1:s 
%
% Solve small system and make v orthogonal to P:
      c = M(k:s,k:s)\f(k:s); 
      v = r - G(:,k:s)*c;
%
% Preconditioning:
      if ( precL )
         if funL
            v = feval( function_M1,v,M1 );
         else
            v = M1\v;
         end
      end
      if ( precU )
         if funU
            v = feval( function_M2,v,M2 );
         else
            v = M2\v;
         end
      end
%
% Compute new U(:,k) and G(:,k), G(:,k) is in space G_j
      U(:,k) = U(:,k:s)*c + om*v;
      if ( funA )
         G(:,k) = feval( function_A,U(:,k),A );
      else
         G(:,k) = A*U(:,k);
      end
%
% Bi-Orthogonalise the new basis vectors: 
      for i = 1:k-1
         alpha =  ( P(:,i)'*G(:,k) )/M(i,i);
         G(:,k) = G(:,k) - alpha*G(:,i);
         U(:,k) = U(:,k) - alpha*U(:,i);
      end
%
% New column of M = P'*G  (first k-1 entries are zero)
      M(k:s,k) = (G(:,k)'*P(:,k:s))';
      if ( M(k,k) == 0 )
         flag = 3;
         return;
      end
%
%  Make r orthogonal to q_i, i = 1..k 
      beta = f(k)/M(k,k);
      r = r - beta*G(:,k);
      x = x + beta*U(:,k);
      normr = norm(r);
      if ( replacement && normr > tolb/mp ) trueres = 1; end;
%
%  Smoothing:
      if ( smoothing )
         t = r_s - r;
         gamma =  (t'*r_s)/(t'*t);
         r_s = r_s - gamma*t;
         x_s = x_s - gamma*(x_s - x);
         normr = norm(r_s);
      end
      resvec = [resvec;normr];
      iter = iter + 1;
      if ( normr < tolb | iter == maxit ) 
         break
      end 
%
% New f = P'*r (first k  components are zero)
      if ( k < s ) 
         f(k+1:s)   = f(k+1:s) - beta*M(k+1:s,k);
      end
   end 
%
   if ( normr < tolb | iter == maxit )
      break
   end
%
% Now we have sufficient vectors in G_j to compute residual in G_j+1
% Note: r is already perpendicular to P so v = r
%
% Preconditioning:
   v = r;
   if ( precL )
      if funL
         v = feval( function_M1,v,M1 );
      else
         v = M1\v;
      end
   end
   if ( precU )
      if funU
         v = feval( function_M2,v,M2 );
      else
         v = M2\v;
      end
   end
%
% Matrix-vector multiplication:
   if ( funA )
      t = feval( function_A,v,A );
   else
      t = A*v;
   end
%
% Computation of a new omega
   om = omega( t, r, angle );
   if ( om == 0 )
      flag = 3;
      return;
   end
%
   r = r - om*t;
   x = x + om*v;
   normr = norm(r);
   if ( replacement && normr > tolb/mp ) trueres = 1; end;
%
%     Residual replacement?
   if ( trueres && normr < normb )
      if funA
         r = b - feval( function_A,x,A );
      else
         r = b - A*x;
      end 
      trueres = 0;
      replacements = replacements+1;
   end
%
%     Smoothing:
   if ( smoothing )
      t    = r_s - r;
      gamma = (t'*r_s)/(t'*t);
      r_s = r_s - gamma*t;
      x_s = x_s - gamma*(x_s - x);
      normr = norm(r_s);
   end
%
   resvec = [resvec;normr];
   iter = iter + 1;

end; %while

if ( smoothing )
   x = x_s;
end

if funA
   relres = norm(b - feval( function_A,x,A ))/normb;
else
   relres = norm(b - A*x)/normb;
end 
if ( relres < tol ) 
   flag = 0;
elseif ( iter == maxit )
   flag = 1;
else
   flag = 2;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function om = omega( t, s, angle )

ns = norm(s);
nt = norm(t);
ts = t'*s;
rho = abs(ts/(nt*ns));
om=ts/(nt*nt);
if ( rho < angle )
   om = om*angle/rho;
end

return
