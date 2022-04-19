% HESTONTESTSPEED - Compare cpu time vs error of various solvers on
%                   the Heston problem.
%
% DESCRIPTION:
%   Plots error against cpu time in the discretized Heston PDE.
%   We only include the ode15s, phipm and Crank--Nicolson methods.
%   For implementation details of the ADI methods the reader is
%   referred to the paper by In't Hout in the references.

% Parameters
K = 100;
Ls = 800;
Lv = 5;
Smax = 200;
Lv = 5;
Vmax = 1;
Ns = 100;
Nv = 50;
Nvpo = Nv+1;

% Set Casee to 0 for the example in In 't Hout (2007), "ADI schemes
% in the numerical solution of the Heston PDE", Numer. Anal. App.
% Math. 936:10-14.
% Set Casee to 1-4 for the examples in In 't Hout and Welfert
% (2009), "Unconditional stability of second-order ADI schemes
% applied to multi-dimensional diffusion equations with mixed
% derivative terms", Appl. Numer. Math. 59(4):677-692.
Casee = 0;

% Heston parameters
if Casee == 0;
  rd = 0.03;
  rf = 0.0;
  nu = 0.2;
  kappa = 2;
  lambda = 0.3;
  rho = 0.8;
  tspan = [0, 1];
elseif Casee == 1;
  rd = 0.025;
  rf = 0.0;
  nu = 0.04;
  kappa = 1.5;
  lambda = 0.3;
  rho = -0.9;
  tspan = [0, 1];
elseif Casee == 2;
  rd = 0.01;
  rf = 0.04;
  nu = 0.12;
  kappa = 3;
  lambda = 0.04;
  rho = 0.6;
  tspan = [0, 1];
elseif Casee == 3;
  rd = 0.03;
  rf = 0.0;
  nu = 0.0707;
  kappa = 0.6067;
  lambda = 0.2928;
  rho = -0.7571;
  tspan = [0, 3];
elseif Casee == 4;
  rd = 0.0507;
  rf = 0.0469;
  nu = 0.06;
  kappa = 2.5;
  lambda = 0.5;
  rho = -0.1;
  tspan = [0, 0.25];
end

% Compute the finite difference discretization 
ds = Ls/(Ns+1);
dv = Lv/(Nv+1);
[L0, L1, L2, L3, L4, L5, b] = FDHeston(ds, dv, Ns, Nv, rd, rf, nu, kappa, lambda, rho);

% The resulting matrix
A = L0+L1+L2+L3+L4+L5;

% -------------------------------------------------------------
% Intial condition
% -------------------------------------------------------------

s = (1:Ns)'*ds;
v = (0:Nv)'*dv;
z = zeros(Ns, 1);
U0 = max(s-K, z);
temp = [];
for j = 1:Nvpo
  temp = [temp; U0];
end;
U0 = temp;

% Set up the correct meshgrid for ploting
I = find(s<=Smax);
s = s(I);
Nps = length(s);
I = find(v<=Vmax);
v = v(I);
Npv = length(v);

% -------------------------------------------------------------
% Exact Solver
% -------------------------------------------------------------

fprintf('------------------------------\n');
fprintf('*** Exact ***\n');
TolMax = 1e-10;
tic;
SolExact = phipm(tspan(end), A, [U0, b], TolMax);
SolExact = reshape(SolExact, Nvpo, Ns);
SolExact = SolExact(1:Nps,1:Npv);
tm = toc;
fprintf('%5.2f cpu time\n', tm);

% Tolerance requirements
Tol = logspace(-1, -6, 7);

for i = 1:length(Tol)

  % -------------------------------------------------------------
  % PHIPM Solver
  % -------------------------------------------------------------

  fprintf('------------------------------\n');
  fprintf('*** PHIPM ***\n');
  tic;
  SolExp = phipm(tspan(end), A, [U0, b], 1e2*Tol(i), false);
  tm = toc;
  SolExp = reshape(SolExp, Nvpo, Ns);
  SolExp = SolExp(1:Nps,1:Npv);
  cpuExp(i) = tm;
  ErrExp(i) = norm(SolExp-SolExact, 2);
  fprintf('%5.2f cpu time\n', tm);
  fprintf('error: %e\n', ErrExp(i));
 
  % -------------------------------------------------------------
  % ode15s Solver
  % -------------------------------------------------------------

  fprintf('------------------------------\n');
  fprintf('*** ode15s ***\n');
  tic;
  odeoptions = odeset('RelTol', 1e-2*Tol(i), 'AbsTol', Tol(i), ...
                      'Jacobian', A);
  Solode15s = ode15s(@rhsheston, tspan, U0, odeoptions, A, b);
  Solode15s = deval(Solode15s, tspan(end));
  tm = toc;
  Solode15s = reshape(Solode15s, Nvpo, Ns);
  Solode15s = Solode15s(1:Nps,1:Npv);
  cpuode15s(i) = tm;
  Errode15s(i) = norm(Solode15s-SolExact, 2);
  fprintf('%5.2f cpu time\n', tm);
  fprintf('error: %e\n', Errode15s(i));
  
end

% -------------------------------------------------------------
% Crank-Nicolson and ADI fixed step solvers
% -------------------------------------------------------------

A0 = L4;
A1 = L0 + L1 + 1/2*L5;
A2 = L2 + L3 + 1/2*L5;
NumSteps = 2.^(8:14);

for i = 1:length(NumSteps)

  % -------------------------------------------------------------
  % Crank-Nicolson Solver
  % -------------------------------------------------------------

  fprintf('------------------------------\n');
  fprintf('*** Crank-Nicolson ***\n');
  tic;
  SolCrank = Crank(tspan, U0, NumSteps(i), A, b);
  tm = toc;
  SolCrank = reshape(SolCrank, Nvpo, Ns);
  SolCrank = SolCrank(1:Nps,1:Npv);
  cpuCrank(i) = tm;
  ErrCrank(i) = norm(SolCrank-SolExact, 2);
  fprintf('%5.2f cpu time\n', tm);
  fprintf('error: %e\n', ErrCrank(i));

  % -------------------------------------------------------------
  % Douglas Solver
  % -------------------------------------------------------------
  
  fprintf('------------------------------\n');
  fprintf('*** Douglas ***\n');
  theta = 1/2;
  tic;
  SolDouglas = Douglas(tspan, U0, NumSteps(i), A0, A1, A2, b, theta);
  tm = toc;
  SolDouglas = reshape(SolDouglas, Nvpo, Ns);
  SolDouglas = SolDouglas(1:Nps,1:Npv);
  cpuDouglas(i) = tm;
  ErrDouglas(i) = norm(SolDouglas-SolExact, 2);
  fprintf('%5.2f cpu time\n', tm);
  fprintf('error: %e\n', ErrDouglas(i));
  
  % -------------------------------------------------------------
  % Hundsdorfer and Verwer Solver
  % -------------------------------------------------------------
  
  fprintf('------------------------------\n');
  fprintf('*** Hundsdorfer and Verwer ***\n');
  theta = 0.3;
  sigma = 1/2;
  tic;
  SolHunds = Hundsdorfer(tspan, U0, NumSteps(i), A0, A1, A2, b, theta, sigma);
  tm = toc;
  SolHunds = reshape(SolHunds, Nvpo, Ns);
  SolHunds = SolHunds(1:Nps,1:Npv);
  cpuHunds(i) = tm;
  ErrHunds(i) = norm(SolHunds-SolExact, 2);
  fprintf('%5.2f cpu time\n', tm);
  fprintf('error: %e\n', ErrHunds(i));

end

figure
clf;
axes('FontSize', 12);
loglog(cpuExp, ErrExp, 'ko-', 'LineWidth', 2);
axis([1e-1 1e2 1e-9 1e-1]);
hold on
loglog(cpuCrank, ErrCrank, 'co-', 'LineWidth', 2);
loglog(cpuDouglas, ErrDouglas, 'go-', 'LineWidth', 2);
loglog(cpuHunds, ErrHunds, 'bo-', 'LineWidth', 2);
loglog(cpuode15s, Errode15s, 'ro-', 'LineWidth', 2);
legend('Phipm', 'Crank-Nicolson', 'Douglas', 'Hundsdorfer', 'Ode15s', ...
       'Location', 'SouthWest');
xlabel('cpu time');
ylabel('maximum error');
