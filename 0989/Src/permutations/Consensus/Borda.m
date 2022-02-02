function[ConsensusRanking] = Borda(DistanceFun,SelPop,NSel)
% [ConsensusRanking] = Borda(DistanceFun,SelPop,NSel)
% Borda:  Calculates the consensus ranking from permutations
%		  finding the average permutation.
%         Note: Only for use with the Mallows. 
%             
% DistanceFun: Distance function to calculate the distance
% SelPop:  Population from which the model is learned
% NSel:    Population size
% OUTPUTS
% ConsensusRanking: The consensus ranking from the population
%
% Created version 12/21/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 02/07/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

AuxConsensusRanking = mean(SelPop); % Calculate the mean from each column

[auxval,min_i] = sort(AuxConsensusRanking); % sort the means, min_i contains the corresponding original index position

[auxval,CR] = sort(min_i); % sort the min_i, CR contains the corresponding original index position 

% The consensus ranking
ConsensusRanking = CR;