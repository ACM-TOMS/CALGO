function[Theta] =  CalculateThetaParameterU_Weighted(ConsensusRanking,Pop,initialTheta,upperTheta,maxit,W)
% [Theta] =  CalculateThetaParameterU_Weighted(ConsensusRanking,Pop,initialTheta,upperTheta,maxit,W)
% CalculateThetaParameterU_Weighted:  Calculates the the Theta Parameter according
% to the given weight W for the Ulam model
%         
% INPUTS
% ConsensusRanking:  Consensus ranking computed for the population
% Pop:  Population
% initialTheta:  Initial theta parameter
% upperTheta: Upper limit for Theta
% maxit: Maximum number of iterations
% W: Weights set for each solution 
%
% OUTPUTS
% Theta: Calculated Theta values.
%
%
% References:
% [1] E. Irurozki, J. Ceberio, B. Calvo. J. A. Lozano. Sampling and learning the Mallows model under the Ulam distance. Technical report EHU-KZAA-TR;2014-04. January 2014.
% [2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011
%
% Created version 20/02/2015. Roberto Santana (roberto.santana@ehu.es)
%
% Last version  20/02/2015. Roberto Santana (roberto.santana@ehu.es)


  [auxval,invConsensus] =sort(ConsensusRanking);
   N = size(Pop,1);
   n = size(Pop,2);
   
composition = Pop(:,invConsensus);
for i=1:N,
 xjs(i,:) = xVector(composition(i,:)); 
end
% To calculate the weighted mean, we first multiply the population with the
% corresponding weight and then divide by the sum of the weights.
xjsmean = sum(bsxfun(@times,xjs,W)) / sum(W); % Calculate the weighted mean from each column
 

Theta = NewtonRaphson(initialTheta,xjsmean,upperTheta,maxit,@CayleyThetaFunction,@CayleyThetaDevFunction,{n});

