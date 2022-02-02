function[kdist] =  kendall_distance(perm1,perm2)
% [kdist] =  kendall_distance(perm1,perm2)
% kendall_distance: Calculate the Kendall-taudistance. Given two permutations, it measures
% the numbers of pairs of items which have opposing ordering.
%    In other words, the minimum number of adjacent transpositions nedded to bring perm1^-1 into perm2^-1.
%
% References:
% [1] P. Diaconis: Group representations in probability and statistics. Institute of Mathematical Statistics Volume11, 1988
% [2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011
%
% Created version 12/10/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 03/25/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es)

[auxval,invperm2] =sort(perm2); %Calculate the inverse of a permutation perm2.
composition = perm1(invperm2); %Calculate the composition of two permutations.
% Sum the adjacent trasnpositions of each position. Sum the descomposition of the composition variable.
% See the vVector definition.
kdist = sum(vVector(composition));


