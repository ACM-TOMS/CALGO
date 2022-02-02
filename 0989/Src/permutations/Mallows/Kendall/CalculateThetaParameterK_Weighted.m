function[Theta] =  CalculateThetaParameterK_Weighted(ConsensusRanking,Pop,initialTheta,upperTheta,maxit,W)
% [Theta] =  CalculateThetaParameterK_Weighted(ConsensusRanking,Pop,initialTheta,upperTheta,maxit,W)
% CalculateThetaParameterK_Weighted:  Calculates the the Theta Parameter according
% to the given weight W for the Kendall model
%         
% INPUTS
% ConsensusRanking:  Consensus ranking computed for the population
% Pop:  Population
% initialTheta:  Initial theta parameter
% upperTheta: Upper limit for Theta
% maxit: Maximum number of iterations
% W: Weights set for each solution 
% OUTPUTS
% Theta: Calculated Theta values.
%
% Created version 06/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 06/05/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

  [auxval,invConsensus] =sort(ConsensusRanking);
   N = size(Pop,1);
   n = size(Pop,2);
   
composition = Pop(:,invConsensus);
for i=1:N,
 vjs(i,:) = vVector(composition(i,:)); 
end
% To calculate the weighted mean, we first multiply the population with the
% corresponding weight and then divide by the sum of the weights.
vjsmean = sum(bsxfun(@times,vjs,W)) / sum(W); % Calculate the weighted mean from each column

 

 Theta = NewtonRaphson(initialTheta,vjsmean,upperTheta,maxit,@KendallThetaFunction,@KendallThetaDevFunction,{n});

