function[ConsensusRanking] = BestPermutation(DistanceFun,SelPop,NSel)
% [ConsensusRanking] = BestPermutation(DistanceFun,SelPop,NSel)
% BestPermutation:  Calculates the consensus ranking from permutations
%		  choosing the permutation with the best fitness (best is the maximum).
%         Note: The permutations are ordered with the best in the first position.
%             
% DistanceFun: Distance function to calculate the distance
% SelPop:  Population from which the model is learned
% NSel:    Population size
% OUTPUTS
% ConsensusRanking: The consensus ranking from the population
%
% Created version 04/08/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/08/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

% Take the first permutation, the first is the best
CR = SelPop(1,:); % CR contains the permutation with the best fitness 

% The consensus ranking
ConsensusRanking = CR;