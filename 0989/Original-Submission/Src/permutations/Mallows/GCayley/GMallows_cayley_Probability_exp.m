function[probability_exp] =  GMallows_cayley_Probability_exp(perm,ConsensusRanking,theta,Psis)
% [probability_exp] = GMallows_cayley_Probability_exp(perm,ConsensusRanking,Psis)
% GMallows_cayley_Probability_exp: Calculates the probability's exponent part of a given permutation in the current model.
%
%
% INPUTS
% perm: a vector that contains the permutation
% ConsensusRanking:  The consensus ranking from the population
% theta: The theta parameter
% Psis: The psi values
% OUTPUTS
% probability_exp: The probability's exponent part of a given permutation in the current model.
%
% References:
% [1] E. Irurozki, B. Calvo, J.A Lozano: Sampling and learning mallows and generalized mallows models under the cayley distance. University of the Basque Country, Tech. Rep., 2013
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: Extending Distance-based Ranking Models in Estimation of Distribution Algorithms. In Proceedings of the Congress on Evolutionary Computation, 2014
%
% Created version 06/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 06/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

% Calculate the probability's exponent part of the permutation in the current model
invConsensus = Invert(ConsensusRanking); %Calculate the inverse of onsensusRanking.
composition = Compose(perm,invConsensus); %Calculate the composition of two permutations.
xjs = xVector(composition);

exponent = sum(theta.*xjs); % Calculate exponent by multipling the psi normalization constant
probability_exp = -exponent;