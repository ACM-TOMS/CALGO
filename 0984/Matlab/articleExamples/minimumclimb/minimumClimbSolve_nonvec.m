function sol = minimumClimbSolve_nonvec(order,K_init,maxK,N_k)
% sol = minimumClimbSolve_nonvec(order,K_init,maxK,N_k)
% solves minimum climb problem using IPOPT and adigator in the vectorized
% mode - derivatives of F(X) computed in non-vectorized mode using the
% black box adigatorGenFiles4Ipopt
% This runs much slower than the vectorized version
% Optimal control problem solved using VERY basic h-method mesh refinement,
% evenly spaced mesh, double number of intervals on each iteration, use
% previous solution as initial guess
%
% Inputs: (all optional, all positive integers)
% order:   first or second order, 1 or 2           (default is 1)
% K_init : inititial mesh number of mesh intervals (default is 4)
% maxK:    maximum number of mesh intervals        (default is 2^10)
% N_K:     number of LGR points per mesh interval  (default is 4)
%
% Notes on time:
% order = 1, all else default, time is approx 25s
% order = 2, all else default, time is approx 175s

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
  % Input structure for adigatorGenFiles4Ipopt
  setup.order = order;
  setup.numvar = nz;
  setup.objective = 'minimumClimbObjective';
  setup.constraint = 'minimumClimbConstraints';
  setup.auxdata = probinfo;
  
  funcs = adigatorGenFiles4Ipopt(setup);
  
  sol(i).tinit = toc(ti);
  
  options.lb = z.min;
  options.ub = z.max;
  options.cl = F.min;
  options.cu = F.max;
  
  options.ipopt.tol = sqrt(eps);
  %options.ipopt.print_level = 0;
  if order == 1
    options.ipopt.hessian_approximation = 'limited-memory';
  end
  
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