function[Theta] =  CalculateThetaParameterGK(ConsensusRanking,Pop,initialTheta,upperTheta,maxit)
% [Theta] =  CalculateThetaParameterGK(ConsensusRanking,Pop,initialTheta,upperTheta,maxit)
% CalculateThetaParameterGK:  Calculates the Theta Parameter for the generalized Kendal model
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
% Created version 04/11/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/15/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

  [auxval,invConsensus] =sort(ConsensusRanking);
   N = size(Pop,1);
   n = size(Pop,2);
   
composition = Pop(:,invConsensus);
for i=1:N,
 vjs(i,:) = vVector(composition(i,:)); 
end
vjsmean = mean(vjs);

 
for j=1:N-1,
 Theta(j) = NewtonRaphson(initialTheta,vjsmean,upperTheta,maxit,@GKendallThetaFunction,@GKendallThetaDevFunction,{n,j});
end
