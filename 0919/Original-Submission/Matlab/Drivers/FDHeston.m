function [L0, L1, L2, L3, L4, L5, b] = FDHeston(ds, dv, Ns, Nv, rd, rf, nu, kappa, lambda, rho)
% FDHESTON - The finite difference discretization of the 2D Heston model.
%
%            u_t = 1/2 s^2 v u_ss + (rd - rf) s u_s + 1/2 lambda^2 v u_vv + 
%                  kappa (nu - v) u_v + rho lambda u_vs - rd s u
%
% DESCRIPTION:
%   Evaluates the linear term of the semi-discretised 2D Heston
%   equation in form u_t = L u + b.
%
% PARAMETERS:
%   ds     - length between grid points in the s domain.
%   dv     - length between grid points in the v domain.
%   Ns     - number of grid points in the s domain.
%   Nv     - number of grid points in the v domain.
%   rd     - domestic interest rate.
%   rf     - foreign interest rate.
%   nu     - mean reversion level.
%   kappa  - rate at which we revert to the mean value.
%   lambda - volatility factor
%   rho    - correlation between underlying Brownian motions.
%
% RETURNS:
%   L0     - discretization of u_ss term.
%   L1     - discretization of u_s term.
%   L2     - discretization of u_v term.
%   L3     - discretization of u_vv term.
%   L4     - discretization of u_sv term.
%   L5     - scaling of u term.
%   b      - boundary effects.

% Save some time
rl = rho*lambda;
hls = 1/2*lambda^2;
r = rd-rf;

% Preallocating space
Nvpo = Nv+1;
Id = speye(Ns*Nvpo, Ns*Nvpo);
L0 = sparse([], [], [], Ns*Nvpo, Ns*Nvpo, 5*Ns*Nvpo);
L1 = sparse([], [], [], Ns*Nvpo, Ns*Nvpo, 5*Ns*Nvpo);
L2 = sparse([], [], [], Ns*Nvpo, Ns*Nvpo, 5*Ns*Nvpo);
L3 = sparse([], [], [], Ns*Nvpo, Ns*Nvpo, 5*Ns*Nvpo);
L4 = sparse([], [], [], Ns*Nvpo, Ns*Nvpo, 5*Ns*Nvpo);
L5 = -rd*Id;

% Used for setting up the main matrix
s = (1:Ns)'*ds;
ss = s.^2;
ss1 = ss;
ss1(end) = 2*ss1(end);
s2 = s;
s2(end) = 0;
v = (0:Nv)'*dv;
z = zeros(Ns, 1);

% Do the finite difference discretization
for j = 1:Nvpo
  L0((j-1)*Ns+1:j*Ns, (j-1)*Ns+1:j*Ns) = ...
      v(j)/(2*ds^2)*spdiags([ss, -2*ss, ss1], -1:1, Ns, Ns)';
  L1((j-1)*Ns+1:j*Ns, (j-1)*Ns+1:j*Ns) = ...
      r/(2*ds)*spdiags([s2, z, -s2], -1:1, Ns, Ns)';
  L3((j-1)*Ns+1:j*Ns, (j-1)*Ns+1:j*Ns) = -2*hls*v(j)/dv^2*speye(Ns, Ns);
end;

for j = 2:Nvpo
  L2((j-1)*Ns+1:j*Ns, (j-2)*Ns+1:(j-1)*Ns) = -kappa*(nu-v(j))/(2*dv)*speye(Ns, Ns);
  L2((j-2)*Ns+1:(j-1)*Ns, (j-1)*Ns+1:j*Ns) = kappa*(nu-v(j-1))/(2*dv)*speye(Ns, Ns);

  L3((j-1)*Ns+1:j*Ns, (j-2)*Ns+1:(j-1)*Ns) = hls*v(j)/dv^2*speye(Ns, Ns);
  L3((j-2)*Ns+1:(j-1)*Ns, (j-1)*Ns+1:j*Ns) = hls*v(j-1)/dv^2*speye(Ns, Ns);
  
  L4((j-1)*Ns+1:j*Ns, (j-2)*Ns+1:(j-1)*Ns) = ...
      rl*v(j)/(4*ds*dv)*spdiags([-s2, z, s2], -1:1, Ns, Ns)';
  L4((j-2)*Ns+1:(j-1)*Ns, (j-1)*Ns+1:j*Ns) = ...
      rl*v(j-1)/(4*ds*dv)*spdiags([s2, z, -s2], -1:1, Ns, Ns)';
end;
L2(1:Ns, 1:Ns) = -3*kappa*nu/(2*dv)*speye(Ns, Ns);
L2(1:Ns, Ns+1:2*Ns) = 4*kappa*nu/(2*dv)*speye(Ns, Ns);
L2(1:Ns, 2*Ns+1:3*Ns) = -kappa*nu/(2*dv)*speye(Ns, Ns);

% The vector b with boundary effects
b0 = s(end)^2/ds*[zeros(Ns-1,1); 1];
temp = [];
for j = 1:Nvpo
  temp = [temp; v(j)*b0];
end;
b0 = temp;
b1 = r*s(end)*[zeros(Ns-1,1); 1];
temp = [];
for j = 1:Nvpo
  temp = [temp; b1];
end;
b1 = temp;
b2 = [zeros(Ns*(Nvpo-1),1); kappa*(nu-v(end))/(2*dv)*s];
b3 = [zeros(Ns*(Nvpo-1),1); hls*v(end)/dv^2*s];
b4 = [zeros(Ns*(Nvpo-1),1); rl*v(end)/(2*dv)*s(1:end-1); 0];
b = b0+b1+b2+b3+b4;
