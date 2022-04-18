function[Theta] =  CalculateThetaParameterU(ConsensusRanking,mindist,ncount,Pop,initialTheta,upperTheta,maxit)
% [Theta] =  CalculateThetaParameterU(ConsensusRanking,Pop,initialTheta,upperTheta,maxit)
% CalculateThetaParameterU:  Computes Theta parameters for Mallow model with Ulam distance 
% 
% INPUTS     
% ConsensusRanking:  Consensus ranking computed for the population
% mindist: 
% ncount: 
% Pop:  Population
% initialTheta:  Initial theta parameter
% upperTheta: Upper limit for Theta
% maxit: Maximum number of iterations
% OUTPUTS
% Theta: Calculated Theta values.
% Created version 04/04/2015. Roberto Santana (roberto.santana@ehu.es) 
%
% Last version 04/20/2015.  Roberto Santana (roberto.santana@ehu.es)  

 N = size(Pop,1);
 n = size(Pop,2);


Theta = NewtonRaphson(initialTheta,mindist,upperTheta,maxit,@UlamThetaFunction,@UlamThetaDevFunction,{n,ncount});

