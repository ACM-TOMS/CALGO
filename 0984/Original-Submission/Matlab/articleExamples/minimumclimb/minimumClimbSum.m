function f = minimumClimbSum(x,lambda,auxdata)
% Get a scalar for second derivative..
F = minimumClimbDynamics(x,auxdata);
f = lambda.'*F(:);