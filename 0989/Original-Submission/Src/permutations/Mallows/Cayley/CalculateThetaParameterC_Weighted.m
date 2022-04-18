function[Theta] =  CalculateThetaParameterC_Weighted(ConsensusRanking,Pop,initialTheta,upperTheta,maxit,W)
% [Psi] =  CalculateThetaParameterC_Weighted(ConsensusRanking,Pop,initialTheta,upperTheta,maxit,W)
% CalculateThetaParameterC_Weighted:  Calculates the the Theta Parameter according
% to the given weight W
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
% Created version 06/05/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 06/05/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

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

