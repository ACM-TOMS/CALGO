function[Theta] =  CalculateThetaParameterK(ConsensusRanking,Pop,initialTheta,upperTheta,maxit)
% [Theta] =  CalculateThetaParameterK(ConsensusRanking,Pop,initialTheta,upperTheta,maxit)
% CalculateThetaParameterK:  Calculates the the Theta Parameter for the Kendal model
% 
% INPUTS     
% ConsensusRanking:  Consensus ranking computed for the population
% Pop:  Population
% initialTheta:  Initial theta parameter
% upperTheta: Upper limit for Theta
% maxit: Maximum number of iterations
% OUTPUTS
% Theta: Calculated Theta values.

%
% Created version 12/23/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

  [auxval,invConsensus] =sort(ConsensusRanking);
   N = size(Pop,1);
   n = size(Pop,2);
   
composition = Pop(:,invConsensus);
for i=1:N,
 vjs(i,:) = vVector(composition(i,:)); 
end
vjsmean = mean(vjs);

 

 Theta = NewtonRaphson(initialTheta,vjsmean,upperTheta,maxit,@KendallThetaFunction,@KendallThetaDevFunction,{n});

