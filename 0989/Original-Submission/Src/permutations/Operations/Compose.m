function[composition] = Compose(perm1,perm2)
% [composition] = Compose(s1,s2)
% Compose:  Calculates the composition of 2 permutation with same size
%             
% perm1:  First permutation
% perm2:  Second permutation
% OUTPUTS
% composition: The composition of 2 permutations
%
% Last version 12/22/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

composition = perm1(perm2);