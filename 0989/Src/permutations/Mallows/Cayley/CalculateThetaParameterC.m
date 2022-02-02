function[Theta] =  CalculateThetaParameterC(ConsensusRanking,Pop,initialTheta,upperTheta,maxit)
% [Psi] =  CalculateThetaParameterC(ConsensusRanking,Pop,initialTheta,upperTheta,maxit)
% CalculateThetaParameterC:  Calculates the the Theta Parameter 
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
% Created version 04/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/04/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

  [auxval,invConsensus] =sort(ConsensusRanking);
   N = size(Pop,1);
   n = size(Pop,2);
   
composition = Pop(:,invConsensus);
for i=1:N,
 xjs(i,:) = xVector(composition(i,:)); 
end
xjsmean = mean(xjs);
 


 Theta = NewtonRaphson(initialTheta,xjsmean,upperTheta,maxit,@CayleyThetaFunction,@CayleyThetaDevFunction,{n});

