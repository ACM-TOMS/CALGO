function[kdist] =  ulam_distance(perm1,perm2)
% [kdist] =  ulam_distance(perm1,perm2)
% ulam_distance: Calculate the Ulam distance. Given two permutations perm1 and perm2, it measures
% the size of the complement of the longest common subsequence between perm1 and perm2
% Alternatively, it can be computed as n minus the length of the longest increasing subsequence (LIS) in
% perm1 * (perm2)^-1.  In this implementation, the alternative way is used to compute the distance
%
% References:
% [1] E. Irurozki, J. Ceberio, B. Calvo. J. A. Lozano. Sampling and learning the Mallows model under the Ulam distance. Technical report EHU-KZAA-TR;2014-04. January 2014.
% [2] J. Ceberio,  A. Mendiburu, J.A Lozano: Introducing the Mallows Model on Estimation of Distribution Algorithms. In Proceedings of International Conference on Neural Information Processing (ICONIP), 2011
%
% Created version 20/02/2015. Roberto Santana (roberto.santana@ehu.es)
%
% Last version  20/02/2015. Roberto Santana (roberto.santana@ehu.es)

[auxval,invperm2] =sort(perm2); %Calculate the inverse of a permutation perm2.
composition = perm1(invperm2); %Calculate the composition of two permutations.

n = size(perm2,2);               
auxSIS = [1,zeros(1,n-1)];  % Auxiliary vector to compute length of the longest increasing subsequence
for i=2:n,
   if(composition(i)>composition(i-1))
     auxSIS(i) = auxSIS(i-1) + 1;
   else
     auxSIS(i) = 1;
   end
end
 
% The distance is the longest increasing subsequence in composition, i.e. the
% maximum value of auxSIS;
kdist = max(auxSIS);


