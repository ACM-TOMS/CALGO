function[fthetadev] =  GKendallThetaDevFunction(theta,newton_params)
% [ftheta] =  GKendallThetaDevFunction(theta,newton_params)
% GKendallThetaDevFunction:  Auxiliary function used by the NewtonRaphson
% optimization algorithm to learn the parameters of the model. In
% Newton-Raphson each new point is computed as  x1=x0-f(x0)/f’(x0) where f'
% is the GKendallThetaDevFunction
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
j = newton_params{2}; %theta index

auxnj = n-j+1; % n-j+1 in the function
fthetadev = (-1)*(exp(theta)/(exp(theta)-1)^2) + ((auxnj^2) * exp(auxnj*theta)) / ((exp(auxnj*theta) -1) ^2);
