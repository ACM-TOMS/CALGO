function[fthetadev] =  CayleyThetaDevFunction(theta,newton_params)
% [ftheta] =  CayleyThetaDevFunction(theta,newton_params)
% CayleyThetaDevFunction:  Auxiliary function used by the NewtonRaphson
% optimization algorithm to learn the parameters of the model. In
% Newton-Raphson each new point is computed as  x1=x0-f(x0)/f’(x0) where f'
% is the CayleyThetaDevFunction
% 
% INPUTS
% theta:  Parameter value (point where the function is evaluated)
% the optimization 
% newton_params: Parameters passed from the optimizer needed to evaluate
% the function
% OUTPUTS
% fthetadev: Calculated Thetadev values.
%
% Created version 03/31/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/04/2014. Josian Santamaria  (jasantamaria003@ikasle.ehu.es)

n = newton_params{1}; %number of variables

vectornjs = [1:n-1]; % 1..n-1
fthetadev = sum( (-vectornjs .* exp(theta)) ./ ((exp(theta) + vectornjs) .^2));
