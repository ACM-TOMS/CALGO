function[VProbs] =  VProbMat(NumbVar,Psis,Theta)
% [VProbs] = VProbMat(NumbVar,Theta)
% VProbMat:  Given a Psi constants vector, it generates the Vj probabilities matrix.
%
% INPUTS
% NumbVar: Number of variables
% Psis: A Psi constants vector
% Theta: The theta parameter, it can ve a single value or a vector according of the following:
%           In the case of Mallows it is a single value
%           In the case of GMallows it is a vector of size NumbVar-1
% OUTPUTS
% VProbs: Generated Vj probabilities matrix.
%
% Created version 06/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 06/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

% In order to make the function compatible with the Mallows and GMallows we create a vector of thetas
% when the Theta value is not a vector (Mallows case, all thetas are the same value) otherwise we copy the Theta vector
if(length(Theta) == 1) % If the Theta variable is a single value
    thetas = Theta * ones(1,NumbVar-1); % Generate the theta vector with NumbVar-1 values of Theta variable
else
    thetas = Theta; % The Theta variable is a vector of size NumbVar-1
end

% Calculate Vj probabilities matrix.
VProbs = zeros(NumbVar-1,NumbVar);
for j=1:NumbVar-1,
    for r=0:NumbVar-j,
        upper=exp((-1) * r * thetas(j));
        lower=Psis(j);
        VProbs(j,r+1) =  upper/lower;
    end
end

