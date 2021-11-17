function f=netestf1(x)
%
% f = netestf1(x)
%
% Helical Valley Function.  This function is primarily for testing
% the nonlinear equations package (see NESOLVE), and is used by NEDEMO.
% The function is taken from Appendix B of "Numerical Methods for
% Unconstrained Optimization and Nonlinear Equations" by Dennis and Schnabel.
%

% n = 3;
f = zeros(3,1);
theta = (1/(2*pi)) * atan(x(2)/x(1));
if (x(1) < 0), theta = theta + .5; end
f(1) = 10 * (x(3) - 10*theta);
f(2) = 10 * (sqrt(x(1)^2 + x(2)^2) - 1);
f(3) = x(3);

