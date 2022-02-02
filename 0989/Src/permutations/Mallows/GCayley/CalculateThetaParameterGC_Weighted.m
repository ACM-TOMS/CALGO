function[Theta] =  CalculateThetaParameterGC_Weighted(ConsensusRanking,Pop,initialTheta,upperTheta,maxit,W)
% [Theta] =  CalculateThetaParameterGC_Weighted(ConsensusRanking,Pop,initialTheta,upperTheta,maxit,W)
% CalculateThetaParameterGC_Weighted:  Calculates the the Theta Parameter according
% to the given weight W for the generalized Mallows model
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
xjsmean(xjsmean==0) = 1/n; % Replace 0 values by 1/n
 
j = (1:n-1); % Range of the thetas (1..problem_size)

 Theta = log(n-j) - log(xjsmean ./ (1-xjsmean));
 
 Theta(Theta < 0) = initialTheta; % Theta must be positive

