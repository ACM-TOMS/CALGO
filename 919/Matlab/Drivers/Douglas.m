function yout = Douglas(tspan, y0, numsteps, A0, A1, A2, b, theta)
% DOUGLAS - ADI method of order 2 for ODEs of the type 
%
%           y_t = A y + b 
%
% SYNOPSIS:
%   YOUT = DOUGLAS(TSPAN, Y0, NUMSTEPS, AO, A1, A2, b, theta) 
%
% PARAMETERS: 
%   TSPAN    - is a vector containing the endpoints of the time integration.
%   Y0       - is a vector of initial conditions.
%   NUMSTEPS - number of steps to be taken, fixed stepsize.
%   A0       - matrix relating to the discretization in S variable.
%   A1       - matrix relating to the discretization in v variable.
%   A2       - matrix relating to the discretization from mixed terms.
%   b        - vector of boundary conditions.
%   theta    - method parameter.
%
% RETURNS:
%   YOUT     - the numerical approximation.

dt = (tspan(2)-tspan(1))/numsteps;
Id = speye(size(A0));
A = A0+A1+A2;

mat1 = Id - theta*dt*A1;
[L1, U1, P1, Q1, R1] = lu(mat1);
mat2 = Id - theta*dt*A2;
[L2, U2, P2, Q2, R2] = lu(mat2);

Up = y0;
for i = 1:numsteps
  Y0 = Up + dt*A*Up + dt*b;
  temp1 = Y0 - theta*dt*A1*Up;
  Y1 = Q1 * (U1 \ (L1 \ (P1 * (R1 \ temp1))));
  temp2 = Y1 - theta*dt*A2*Up;
  Up = Q2 * (U2 \ (L2 \ (P2 * (R2 \ temp2))));
end 
yout = Up;
