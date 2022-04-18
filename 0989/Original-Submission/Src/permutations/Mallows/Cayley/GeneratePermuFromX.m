function [permu] = GeneratePermuFromX(x,NumbVar)
% [permu] = GeneratePermuFromX(x,NumbVar)
% GeneratePermuFromX:  Generate a random permutation from a given x vector, that 
% contains the decomposition of the permutation.
%
% INPUTS
% x: a vector that contains,the decomposition of the Cayley distance of a permutation with the permutation identity.
%    In other words, the decomposition of the permutation.
% NumbVar:   Number of variables
% OUTPUTS
% permu: Generated permutation
%
% References:
% [1] E. Irurozki, B. Calvo, J.A Lozano: Sampling and learning mallows and generalized mallows models under the cayley distance. University of the Basque Country, Tech. Rep., 2013
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: Extending Distance-based Ranking Models in Estimation of Distribution Algorithms. In Proceedings of the Congress on Evolutionary Computation, 2014
%
% Created version 03/27/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 06/12/2016. Roberto Santana (roberto.santana@ehu.es) 

permu = [1:NumbVar];   % Init permutation as [1 2 3 4 … n]
for pos=1:NumbVar-1,
   if(x(pos) == 1)  
     random = pos+1+randi([0,NumbVar-pos-1],1);    % random \in (pos+1)+[0..NumbVar-pos-1] = [pos+1..NumbVar]
     aux = permu(random); 		                   % swap positions ‘random’ and ‘pos’ in permu          
     permu(random) = permu(pos);
     permu(pos) = aux;
   end    
end






