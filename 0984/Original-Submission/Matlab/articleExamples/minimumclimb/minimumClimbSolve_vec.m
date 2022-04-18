function [sol,tinit] = minimumClimbSolve_vec(order,K_init,maxK,N_k)
% [sol, tinit] = minimumClimbSolve_vec(order,K_init,maxK,N_k)
% solves minimum climb problem using IPOPT and adigator in the vectorized
% mode - derivatives of F(X) computed in vectorized mode, derivatives of
% NLP built by hand code from those of F(X)
% Optimal control problem solved using VERY basic h-method mesh refinement,
% evenly spaced mesh, double number of intervals on each iteration, use
% previous solution as initial guess
%
% Inputs: (all optional, all positive integers)
% order:   first or second order, 1 or 2           (default is 2)
% K_init : inititial mesh number of mesh intervals (default is 4)
% maxK:    maximum number of mesh intervals        (default is 2^10)
% N_K:     number of LGR points per mesh interval  (default is 4)
%
% Outputs:
% sol:     solution structure
% tinit:   initial time spent

if ~exist('N_k','var')
  N_k = 4;
end
if ~exist('maxK','var')
  maxK = 2^10;
end
if ~exist('K_init','var')
  K_init = 4;
end
if ~exist('order','var')
  order = 2;
end

% solve the minimum climb problem with IPOPT and ADiGator
ipoptexist =  exist('ipopt','file')==3;
if ~ipoptexist
addpath ipopt
end

% get auxdata
auxdata = minimumClimbAuxdata();

% Y = state, U = control, tf = final time
% z = [Y(:), U(:), YNp1(:), tf]
ny = 3; nu = 1; nx = ny+nu;

ti = tic;
% Differentiate Minimum Climb Dynamics
ax = adigatorCreateDerivInput([Inf nx],'X'); % 3 states 1 control
opts = adigatorOptions('overwrite',1,'comments',0);
axdot_X = adigator('minimumClimbDynamics',{ax,auxdata},'minimumClimbDynamics_X',opts);
ifx = axdot_X{1}.deriv.nzlocs(:,1);
jfx = axdot_X{1}.deriv.nzlocs(:,2);
if order == 2
  ax2.f = ax; ax2.dX = adigatorCreateAuxInput([Inf nx]);
  axdot_XX = adigator('minimumClimbDynamics_X',{ax2,auxdata},'minimumClimbDynamics_XX',opts);
  ifxx = ifx(axdot_XX{1}.dX.deriv.nzlocs(:,1));
  jfxx = jfx(axdot_XX{1}.dX.deriv.nzlocs(:,1));
  kfxx = axdot_XX{1}.dX.deriv.nzlocs(:,2);
end
tinit = toc(ti);

sol = struct();
KK = 2.^(floor(log2(K_init)):floor(log2(maxK)));
i  = 0;
for K = KK
  i = i+1;
  if i == 1
    % Setup Problem on Initial Mesh
    [z,F,probinfo] = setupMinimumClimb(N_k,K);
  else
    % Use previous solution as initial guess for the rest..
    [z,F,probinfo] = setupMinimumClimb(N_k,K,sol(i-1));
  end
  nz = length(z.guess);
  N = probinfo.LGR.nLGR;
  probinfo.auxdata = auxdata;
  ti = tic;
  % -------------- Jacobian Projections/Sparsity Patterns
  % Project vectorized locations
  [iFX,jFX] = adigatorProjectVectLocs(N,ifx,jfx);
  probinfo.map.iFX = iFX;
  probinfo.map.jFX = jFX;
  
  % Get d(D*[Y;YNp1])/dz
  jDv = [probinfo.map.state; probinfo.map.stateNp1];
  dVdz = sparse(1:(N+1)*ny,jDv(:),ones(ny*(N+1),1),ny*(N+1),nz);
  dDYdz = reshape(probinfo.LGR.Dmatrix*reshape(dVdz,N+1,ny*nz),N*ny,nz);
  
  probinfo.dDYdz = dDYdz;
  
  dFdz = [sparse(iFX,jFX,ones(size(iFX)),N*ny,nz-1),sparse(1:N*ny,ones(1,N*ny),ones(1,N*ny),N*ny,1)];
  Jpat = dDYdz+dFdz;
  
  if order == 2
    % ------------ Hessian Projections/Sparsity Patterns
    
    [iFXX,jFXX,kFXX] = adigatorProjectVectLocs(N,ifxx,jfxx,kfxx);
    % F_XX is N*ny by N*nx N*nx - roll second two to multiply through by lambda
    jkFXX = sub2ind([N*nx, N*nx],jFXX,kFXX);
    
    
    dLdXdX = sparse(jFXX,kFXX,ones(size(kFXX)),N*nx,N*nx);
    dLdXdX = spones(dLdXdX);
    
    % want to use compression when finding d^2(lambda.'*F)/dx^2..
    dFdXdX = sparse(iFXX,jkFXX,1:numel(iFXX),N*ny,N*nx*N*nx);
    dFdXdXc = dFdXdX(:,unique(jkFXX));
    [iFXXc,jFXXc,kFXXc] = find(dFdXdXc);
    % need to re-order..
    [~,kFXXc_order] = sort(kFXXc);
    iFXXc = iFXXc(kFXXc_order);
    jFXXc = jFXXc(kFXXc_order);

    
    probinfo.map.iFXXc  = iFXXc;
    probinfo.map.jFXXc  = jFXXc;
    probinfo.map.nFXXc  = size(dFdXdXc,2);
    [iLXX,jLXX] = find(dLdXdX);
    probinfo.map.iLXX = iLXX;
    probinfo.map.jLXX = jLXX;
    
    dLdYNp1dz = sparse([],[],[],ny,nz);
    dLdzdt = [sum(dFdz(:,1:end-1),1),0];
    Hpat = tril([[dLdXdX,sparse([],[],[],N*nx,nz-N*nx)];dLdYNp1dz;dLdzdt]);
  end
  sol(i).tinit = toc(ti);
  % ----------------- IPOPT setup ---------------
  
  funcs.objective = @(z)minimumClimbObjective(z,probinfo);
  funcs.gradient  = @grdwrap;
  funcs.constraints = @(z)minimumClimbConstraints(z,probinfo);
  funcs.jacobian = @(z)jacwrap(z,probinfo);
  funcs.jacobianstructure = @()Jpat;
  if order == 2
    funcs.hessian = @(z,sigma,lambda)heswrap(z,sigma,lambda,probinfo);
    funcs.hessianstructure = @()Hpat;
  else
     options.ipopt.hessian_approximation = 'limited-memory';
  end
  
  options.lb = z.min;
  options.ub = z.max;
  options.cl = F.min;
  options.cu = F.max;
  
  options.ipopt.tol = sqrt(eps);
  %options.ipopt.print_level = 0;
  
  % Call IPOPT
  ti = tic;
  [z1,info] = ipopt(z.guess,funcs,options);
  
  % Extract solution, some other stuff
  sol(i).solvetime  = toc(ti);
  sol(i).state = z1(probinfo.map.state);
  sol(i).control = z1(probinfo.map.control);
  sol(i).stateNp1 = z1(probinfo.map.stateNp1).';
  sol(i).tf = z1(probinfo.map.tf);
  sol(i).tau = probinfo.LGR.tau;
  sol(i).NLP.njac = info.eval.jacobian;
  sol(i).NLP.nhes = info.eval.hessian;
  sol(i).nLGR = N;
end

if ~ipoptexist
rmpath ipopt
end

end

function g = grdwrap(z)
% Gradient wrapper
nz = length(z);
g = [zeros(1,nz-1), 1];
end

function J = jacwrap(z,probinfo)
% Jacobian wrapper

Y = z(probinfo.map.state); 
U = z(probinfo.map.control);
tf = z(probinfo.map.tf);

N = probinfo.LGR.nLGR;
nx = 4;
ny = 3;

% Build X = [Y U]
X = [Y, U];
auxdata = probinfo.auxdata;

% Get Dynamics
X_X.f = X;
X_X.dX = ones(size(X));
F = minimumClimbDynamics_X(X_X,auxdata);

% z = [Y(:); U(:); YNp1(:); tf]
% Defect Constraints
%C = D*[Y;YNp1] - (tf/2)*F;
% let A = D*[Y;YNp1], B = (tf/2)*F;
% C = A-B;

% Linear portion of derivative
dAdz = probinfo.dDYdz;

% Nonlinear portion
%dBdX = tf/2*dFdX
dBdX = sparse(probinfo.map.iFX,probinfo.map.jFX,(tf/2)*F.dX,N*ny,N*nx);
%dBdYNp1 = 0
dBdYNp1 = sparse([],[],[],N*ny,ny);
%dBdtf = F/2
dBdtf = F.f(:)/2;

J = dAdz - [dBdX dBdYNp1 dBdtf];
end

function H = heswrap(z,sigma,lambda,probinfo)
% Hessian wrapper

Y = z(probinfo.map.state); 
U = z(probinfo.map.control);
tf = z(probinfo.map.tf);

N = probinfo.LGR.nLGR;
nx = 4;
ny = 3;

nz = length(z);

% Build X = [Y U]
X = [Y, U];
auxdata = probinfo.auxdata;

% Get Dynamics
X_X.f = X;
X_X.dX = ones(size(X));
F = minimumClimbDynamics_XX(X_X,auxdata);

% z = [Y(:); U(:); YNp1(:); tf]
% Defect Constraints
% C = D*[Y;YNp1] - (tf/2)*F;
% B = tf/2*F


% Nonlinear portion
%dBdX = tf/2*dFdX
FXX = sparse(probinfo.map.iFXXc,probinfo.map.jFXXc,-(tf/2)*F.dXdX,N*ny,probinfo.map.nFXXc);
lambda_FXX = lambda.'*FXX;
dLdXdX = sparse(probinfo.map.iLXX,probinfo.map.jLXX,lambda_FXX,N*nx,N*nx);
dLdXdX = 1/2*(dLdXdX+dLdXdX.'); % average
% only need to worry about lower triangular for ipopt..


dLdXdtf = lambda.'*sparse(probinfo.map.iFX,probinfo.map.jFX,-F.dX./2,N*ny,N*nx);

H = [[tril(dLdXdX), sparse([],[],[],N*nx,nz-(N*nx))];sparse([],[],[],ny,nz);[dLdXdtf, sparse([],[],[],1,ny+1)]];
end