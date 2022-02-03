function[ConsensusRanking,min_dist] = setMedianPermutationUlam(DistanceFun,SelPop,NSel)
% [ConsensusRanking] = setMedianPermutationUlam(DistanceFun,SelPop,NSel)
% setMedianPermutationUlam:  Calculates the consensus ranking from permutations
%		  choosing the solution that minimizes the distance to the sample. 
%             
% INPUTS
% DistanceFun: Distance function to calculate the distance
% SelPop:  Population from which the model is learned
% NSel:    Population size
% OUTPUTS
% ConsensusRanking: The consensus ranking from the population
%
% Created version 01/29/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 02/07/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

%Calculate the distances between all pairs of permutations
distances = zeros(NSel);
for i=1:NSel,
    for j=i+1:NSel,
        dist = eval([DistanceFun,'(SelPop(i,:),SelPop(j,:))']);
        distances(i,j) = dist;
    end
end

distances = distances + distances';

AuxConsensusRanking = sum(distances); % Calculate the sum from each column


[value, index]=min(AuxConsensusRanking); % calculate the minimum distance, index contains the position of the minimum value

CR = SelPop(index,:); % CR contains the permutation that minimizes the distance

% The consensus ranking
ConsensusRanking = CR;
min_dist=value/NSel;
