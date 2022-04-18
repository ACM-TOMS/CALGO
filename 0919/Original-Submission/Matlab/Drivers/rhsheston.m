function f = rhsheston(t, y, A, b);
% RHSHESTON - Computes the full rhs of the discretised Heston PDE
%
%             y_t = A y + b 
%
% DESCRIPTION:
%   This is used for the matlab function ode15s.
%
% PARAMETERS: 
%   t    - is a scalar of current time.
%   y    - is a vector of solution values.
%   A    - constant matrix.
%   b    - vector of boundary conditions.
%
% RETURNS:
%   f    - the rhs.

f = A*y+b;
