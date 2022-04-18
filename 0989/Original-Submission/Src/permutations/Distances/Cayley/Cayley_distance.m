function[cdist] =  Cayley_distance(perm1,perm2)
% [cdist] =  Cayley_distance(perm1,perm2)
% Cayley_distance: Calculate the Cayley distance. Given two permutations, it measures
% the numbers of swaps (not necessary adjacent) nedded to convert perm1 into perm2.
%
% References:
% [1] E. Irurozki, B. Calvo, J.A Lozano: Sampling and learning mallows and generalized mallows models under the cayley distance. University of the Basque Country, Tech. Rep., 2013
% [2] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: Extending Distance-based Ranking Models in Estimation of Distribution Algorithms. In Proceedings of the Congress on Evolutionary Computation, 2014
%
% Created version 04/01/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 04/01/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es)

[auxval,invperm2] =sort(perm2); %Calculate the inverse of a permutation perm2.
composition = perm1(invperm2); %Calculate the composition of two permutations.
% Sum the swaps (not necessary adjacent) of each position. Sum the descomposition of the composition variable.
% See the xVector definition.
cdist = sum(xVector(composition));


