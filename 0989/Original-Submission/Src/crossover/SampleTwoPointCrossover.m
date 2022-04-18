function [NewPop] = SampleTwoPointCrossover(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% [NewPop] = SampleTwoPointCrossover(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params)
% SampleTwoPointCrossover: Applies the two-point crossover recombination as
%                          specified by the model and the calls the mutation operator 
%               
% INPUTS
% NumbVar:   Number of variables
% model:      model{1} contains a matrix with one row for each individual, and two columns.
%             The first column is the list of  the first parents. The second column, the list of second parents. 
%             They may coincide, but this is highly unlikely and will produce not harm
%             model{2} contains,  for each pair of solutions, the point
%             where crossover is applied to.
% Card:      Vector with the dimension of all the variables. 
% AuxPop:    Auxiliary (selected) population that it is used here to do
%            crossover between solutions
% AuxFunVal: Evaluation of the data set (required for some sampling algorithms, not for this one)
% sampling_params{1}(1) = N: Number of generated individuals 
% OUTPUTS
% NewPop: Sampled individuals
%
% Last version 11/7/2013. Roberto Santana (roberto.santana@ehu.es)  


[NewPop] = CXTwoPoint(NumbVar,model,Card,AuxPop,AuxFunVal,sampling_params);

mutOperator = char(cellstr(sampling_params{1}(2)));      % Type of mutation to be applied 

for i=3:size(sampling_params{1},2),                  % We assume that all parameters from 2 to the end are mutation parameters and pass them to mutation
  mutation_params{1}(i-2) = sampling_params{1}(i);
end

[NewPop]  = eval([mutOperator,'(NumbVar,Card,NewPop,mutation_params)']);



