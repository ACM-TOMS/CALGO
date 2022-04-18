function[XProbs] =  XProbMat(NumbVar,Psis,Theta)
% [XProbs] = XProbMat(NumbVar,Psis,Theta)
% XProbMat:  Given a Psi constants vector it generates the Xj probabilities matrix.
%
% INPUTS
% NumbVar: Number of variables
% Psis: A Psi constants vector
% Theta: The theta parameter, it can ve a single value or a vector according of the following:
%           In the case of Mallows it is a single value
%           In the case of GMallows it is a vector of size NumbVar-1
% OUTPUTS
% XProbs: Generated Xj probabilities matrix.
%
% Created version 06/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 06/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Calculate Xj probabilities matrix.
XProbs = zeros(NumbVar-1,1);
for j=1:NumbVar-1,
    XProbs(j) = 1/Psis(j);
end

