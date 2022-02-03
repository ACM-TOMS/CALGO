function [NewPop] = SampleTransposition(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = SampleTransposition(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% SampleTransposition: 
%                              Applies a transposition operator to a population of individuals. 
%                              A transposition consists of randomly selecting a position i in the individual, an
%                              offset a, and a length l. The individual is then modified by setting 
%                              x_{i+j} = x_{i+a-j} for j={0,...,l-1}.The parameters for the transposition are received  a model 
%                              constructed with function LearnSymmetricTransposition.m
%                              After applying the transposition the mutation operator contained in
%                              sampling_params{1}(1) is called.
% INPUTS
% NumbVar:   Number of variables
% model:     Contains four parts specifying how the transposition operator
%            is applied.  model{1} stores the indexes of individuals where
%            transposition operator will be applied. model{2} the  length of each
%            transposition. model{3} the initial location of each transposition.  
%            model{4} the offSets of each transposition
% Card:      Vector with the dimension of all the variables. 
% AuxPop:    Auxiliary (selected) population that it is used here to do
%            crossover between solutions
% AuxFunVal: Evaluation of the data set (required for some sampling algorithms, not for this one)
% sampling_params{1}(1) = N: Number of generated individuals 
% OUTPUTS
% NewPop: Sampled individuals
%
% 
% Last version 11/7/2013. Roberto Santana (roberto.santana@ehu.es) following 
% Symmetry in evolutionary and estimation of distribution algorithms.
% (Santana, McKay, Lozano: 2013)  2013 IEEE Congress on Evolutionary
% Computation (CEC).  http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6557811&tag=1


[NewPop] = CXTransposition(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params);

mutOperator = char(cellstr(sampling_params{1}(2)));      % Type of mutation to be applied 

for i=3:size(sampling_params{1},2),                  % We assume that all parameters from 2 to the end are mutation parameters and pass them to mutation
  mutation_params{1}(i-2) = sampling_params{1}(i);
end

[NewPop]  = eval([mutOperator,'(NumbVar,Card,NewPop,mutation_params)']);



 