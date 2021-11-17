function [vals] = Eval_uBQP(ind)
% [vals] = Eval_uBQP(ind)(ind)
%  Eval_uBQP(ind): Evaluates one configuration of the multi_objective uBQP problem for 
%                  Requires that the structure for each objective has been
%                  defined as global structure. It can be read with the LodaUBPInstance.m
%                  function 
% INPUT
% ind: the individual (vector) which represents a matrix of values for the spins
% all_obj: Structures with the interactions between variables pairs 
% OUTPUT
% vals: The evaluation of the uBQP function in each objective for the individual
%
% Last version 6/09/2015. Roberto Santana (roberto.santana@ehu.es)

global all_obj

n_obj = size(all_obj,2);
vals = 0;

for l=1:n_obj,
  vals(l) = 0;
  for i=1:size(all_obj{l},1)
    vals(l) = vals(l) + ind(all_obj{l}(i,1))*ind(all_obj{l}(i,2))*all_obj{l}(i,3);
  end
end  
   





