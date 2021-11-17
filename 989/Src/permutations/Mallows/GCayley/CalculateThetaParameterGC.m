function[Theta] =  CalculateThetaParameterGC(ConsensusRanking,Pop,initialTheta,upperTheta,maxit)
% [Theta] =  CalculateThetaParameterGC(ConsensusRanking,Pop,initialTheta,upperTheta,maxit)
% CalculateThetaParameterGC:  Calculates the the Theta Parameter for the generalized Cayley model
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
% Created version 04/15/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/15/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

  [auxval,invConsensus] =sort(ConsensusRanking);
   N = size(Pop,1);
   n = size(Pop,2);
   
composition = Pop(:,invConsensus);
for i=1:N,
 xjs(i,:) = xVector(composition(i,:)); 
end
xjsmean = mean(xjs);
xjsmean(xjsmean==0) = 1/n; % Replace 0 values by 1/n
 
j = (1:n-1); % Range of the thetas (1..problem_size)

 Theta = log(n-j) - log(xjsmean ./ (1-xjsmean));
 
 Theta(Theta < initialTheta) = initialTheta; % Theta must be positive

