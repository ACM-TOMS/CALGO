function [perm] = GeneratePermuAtLIS(NumbVar,shape) 
%  [perm] = GeneratePermuAtLIS(NumbVar,shape)
% GeneratePermuAtLIS: Generates a permutation at a given distance using a
% Ferrer shape 
% A Ferrers diagram with shape lambda is a set of cells as shown  where row
% i has lambda_i cells.
%         
% INPUTS
% NumbVar: Number of variables
% shape: 
%
% OUTPUTS
% perm: Generated permutation
%
%
% Created version 04/04/2015. Roberto Santana (roberto.santana@ehu.es) 
%
% Last version 04/20/2015.  Roberto Santana (roberto.santana@ehu.es)  
%
% References:
% [1] 
% [2] 

 [YoungT1] =  RandomTableau(NumbVar,shape);
 [YoungT2] =  RandomTableau(NumbVar,shape);
 
 [perm] =  GenPermFromYoungTableaux(YoungT1,YoungT2);

 

 