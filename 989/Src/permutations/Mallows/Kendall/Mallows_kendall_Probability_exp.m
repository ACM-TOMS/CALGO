function[probability_exp] =  Mallows_kendall_Probability_exp(perm,ConsensusRanking,theta,Psis)
% [probability_exp] = Mallows_kendall_Probability_exp(perm,ConsensusRanking,Psis)
% Mallows_kendall_Probability_exp:  Calculates the probability's exponent part of a given permutation 
% in the current model.
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
% [1] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: A Distance-based Ranking Model Estimation of Distribution Algorithm for the Flowshop Scheduling Problem. IEEE Transactions on Evolutionary Computation, 2013
%
% Created version 05/23/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 05/23/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

% Calculate the probability's exponent part of the permutation in the current model
dist = kendall_distance(perm,ConsensusRanking); % Calculate the distance

probability_exp = -theta*dist;