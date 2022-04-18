function[probability] =  GMallows_kendall_Probability(perm,ConsensusRanking,theta,Psis)
% [probability] = GMallows_kendall_Probability(perm,ConsensusRanking,Psis)
% GMallows_kendall_Probability:  Calculates the probability of a given permutation in the current model.
%
% INPUTS
% perm: a vector that contains the permutation
% ConsensusRanking:  The consensus ranking from the population
% theta: The theta parameter
% Psis: The psi values
% OUTPUTS
% probability: The probability of a given permutation in the current model.
%
% References:
% [1] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: A Distance-based Ranking Model Estimation of Distribution Algorithm for the Flowshop Scheduling Problem. IEEE Transactions on Evolutionary Computation, 2013
%
% Created version 06/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 06/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

% Calculate the probability of the permutation in the current model
invConsensus = Invert(ConsensusRanking); %Calculate the inverse of onsensusRanking.
composition = Compose(perm,invConsensus); %Calculate the composition of two permutations.
vjs = vVector(composition);

psi_prod = prod(Psis); % Calculate the psi normalization constant
exponent = sum(theta.*vjs); % Calculate exponent by multipling the psi normalization constant
probability = exp(-exponent)/psi_prod;