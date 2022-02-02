function yout = Crank(tspan, y0, numsteps, A, b)
% CRANK - Crank-Nicolson method of order 2 for ODEs of the type 
%
%          y_t = A y + b 
%
% PARAMETERS: 
%   tspan    - is a vector containing the endpoints of the time integration.
%   y0       - is a vector of initial conditions.
%   numsteps - number of steps to be taken, fixed stepsize.
%   A        - constant matrix.
%   b        - vector of boundary conditions.
%
% RETURNS:
%   yout     - the numerical approximation.

dt = (tspan(2)-tspan(1))/numsteps;
Id = speye(size(A));

% Comment this in when using recent matlab
mat = Id - dt*A/2;
[L, U, P, Q, R] = lu(mat);

Up = y0;
for i = 1:numsteps
  temp = Up + dt*A*Up/2 + dt*b;
  Up = Q * (U \ (L \ (P * (R \ temp))));
end 
yout = Up;
