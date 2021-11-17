function[probability] =  Mallows_Ulam_Probability(perm,ConsensusRanking,theta,VProbs)
% [probability] = Mallows_Ulam_Probability(perm,ConsensusRanking,Psis)
% Mallows_Ulam_Probability:  Calculates the probability of a given permutation in the current model.
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
% [1] E. Irurozki, B. Calvo, J.A Lozano: Sampling and learning mallows and generalized mallows models under the Ulam distance. University of the Basque Country, Tech. Rep., 2013
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: Extending Distance-based Ranking Models in Estimation of Distribution Algorithms. In Proceedings of the Congress on Evolutionary Computation, 2014
%
%
% Created version 02/12/2015. Roberto Santana (roberto.santana@ehu.es) 
%
% Last version 04/20/2015. Roberto Santana (roberto.santana@ehu.es)  


% Calculate the probability of the permutation in the current model
dist = ulam_distance(perm,ConsensusRanking); % Calculate the distance
probability = exp(-dist * theta)/ VProbs(end);
  
  if probability==0,
    printf('Numerical error in Ulam-Probability!!');
    exit()
  end


