function[newvector] =  vVector(vector)
% [newvector] = vVector(vector)
% vVector:  Given a permutation, it generates the decomposition of the permutation.
% The decomposition of the Kendall-tau distance of a permutation with the permutation identity.
%    In other words, the decomposition of the permutation.
%
% INPUTS
% vector: a vector that contains the permutation
%
% OUTPUTS
% newvector: Generated decomposition of the permutation.
%
% References:
% [1] J. Ceberio, E. Irurozki, A. Mendiburu, J.A Lozano: A Distance-based Ranking Model Estimation of Distribution Algorithm for the Flowshop Scheduling Problem. IEEE Transactions on Evolutionary Computation, 2013
%
% Created version 12/15/2013. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 
%
% Last version 03/25/2014. Josian Santamaria (jasantamaria003@ikasle.ehu.es) 

% Calculate the decomposition of the permutation, calculating the number of positions of the permutation
% on the right of i that have values smaller than vector(i)
for i=1:length(vector)-1,
    newvector(i) = sum(vector(i+1:end)<vector(i));
end 