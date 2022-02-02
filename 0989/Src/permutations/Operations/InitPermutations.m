function [NewPop] = InitPermutations(NumbVar,PopSize,Card,seeding_params)
% [NewPop] = InitPermutations(NumbVar,PopSize,Card,sampling_params)
% InitPermutations:         Random Initialization of a population
%                           comprising permutations
% INPUTS
% NumbVar: Number of variables
% Card: For discrete variables:    Vector with the dimension of all the variables. 
% OUTPUTS
% NewPop: Sampled individuals
%

for i=1:PopSize,                      % Each individual is a random permutation
   NewPop(i,:) =  randperm(NumbVar);
end
 
